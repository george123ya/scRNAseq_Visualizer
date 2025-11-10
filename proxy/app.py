from flask import Flask, request, jsonify
import fsspec
import zarr
import json
import os
import pandas as pd
import pyarrow.parquet as pq
import pyarrow as pa
from datetime import datetime
import numpy as np
import logging
import tempfile
import shutil

app = Flask(__name__)
logging.basicConfig(level=logging.DEBUG)

# Load credentials
def load_credentials():
    cred_file = os.path.join(os.path.dirname(__file__), 'credentials.json')
    if not os.path.exists(cred_file):
        raise FileNotFoundError("credentials.json not found.")
    with open(cred_file, 'r') as f:
        return json.load(f)

CREDENTIALS = load_credentials()

def get_filesystem(credential_id='default'):
    creds = CREDENTIALS.get(credential_id)
    if not creds:
        raise ValueError(f"Credential ID {credential_id} not found")
    return fsspec.filesystem('s3', 
                           key=creds.get('aws_access_key_id'), 
                           secret=creds.get('aws_secret_access_key'))

def generate_version_prefix(user_id, custom_prefix=None):
    if custom_prefix:
        return custom_prefix
    timestamp = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    return f"{user_id}_{timestamp}"

def delete_version(fs, zarr_url, version_prefix):
    version_path = f"{zarr_url}/versions/{version_prefix}"
    if fs.exists(version_path):
        app.logger.debug(f"Deleting existing version: {version_path}")
        fs.rm(version_path, recursive=True)

def save_version_metadata(fs, zarr_url, version_prefix, metadata):
    zattrs_path = f"{zarr_url}/versions/{version_prefix}/.zattrs"
    app.logger.debug(f"Saving metadata to: {zattrs_path}")
    with fs.open(zattrs_path, 'w') as f:
        json.dump(metadata, f, indent=2)

def save_parquet_data(fs, zarr_url, version_prefix, data_dict, subfolder="results"):
    for key, data in data_dict.items():
        parquet_path = f"{zarr_url}/versions/{version_prefix}/{subfolder}/{key}.parquet"
        dir_path = f"{zarr_url}/versions/{version_prefix}/{subfolder}"
        if not fs.exists(dir_path):
            fs.mkdir(dir_path)
        
        if isinstance(data, list):
            try:
                df = pd.DataFrame(data)
            except Exception as e:
                raise ValueError(f"Failed to convert {key} to DataFrame: {str(e)}")
        else:
            df = data
        
        app.logger.debug(f"Saving Parquet: {parquet_path}")
        table = pa.Table.from_pandas(df)
        with fs.open(parquet_path, 'wb') as f:
            pq.write_table(table, f)

def save_zarr_array(fs, zarr_url, version_prefix, array_name, data, subfolder="masks"):
    zarr_path = f"{zarr_url}/versions/{version_prefix}/{subfolder}/{array_name}"
    dir_path = f"{zarr_url}/versions/{version_prefix}/{subfolder}"
    if not fs.exists(dir_path):
        fs.mkdir(dir_path)
    
    data_array = np.array(data) if isinstance(data, list) else data
    app.logger.debug(f"Saving Zarr array: {zarr_path}")
    store = fs.get_mapper(zarr_path)
    z = zarr.create(shape=data_array.shape, chunks=data_array.shape, dtype=data_array.dtype, store=store)
    z[:] = data_array

def update_global_index(fs, zarr_url, version_prefix, user_id, metadata):
    bucket = zarr_url.split("://")[1].split('/')[0]
    index_path = f"s3://{bucket}/versions/index.parquet"
    timestamp = datetime.utcnow().isoformat()
    new_entry = pd.DataFrame([{
        'zarr_url': zarr_url,
        'version_prefix': version_prefix,
        'user_id': user_id,
        'timestamp': timestamp,
        'metadata_keys': list(metadata.keys())
    }])
    
    app.logger.debug(f"Updating global index: {index_path}")
    table = pa.Table.from_pandas(new_entry)
    if fs.exists(index_path):
        existing = pq.read_table(index_path, filesystem=fs)
        updated = pa.concat_tables([existing, table])
        with fs.open(index_path, 'wb') as f:
            pq.write_table(updated, f)
    else:
        with fs.open(index_path, 'wb') as f:
            pq.write_table(table, f)

@app.route('/add_version', methods=['POST'])
def add_version():
    try:
        data = request.json
        zarr_url = data['zarr_url']
        user_id = data.get('user_id', 'anonymous')
        version_prefix = data.get('version_prefix')
        metadata = data.get('metadata', {})
        parquet_data = data.get('parquet_data', {})
        zarr_arrays = data.get('zarr_arrays', {})
        credential_id = data.get('credential_id', 'default')

        version_prefix = generate_version_prefix(user_id, version_prefix)
        if not zarr_url:
            return jsonify({'error': 'Missing zarr_url'}), 400

        fs = get_filesystem(credential_id)
        bucket = zarr_url.split("://")[1].split('/')[0]
        versions_dir = f"{zarr_url}/versions"
        if not fs.exists(versions_dir):
            fs.mkdir(versions_dir)

        delete_version(fs, zarr_url, version_prefix)
        save_version_metadata(fs, zarr_url, version_prefix, metadata)
        if parquet_data:
            save_parquet_data(fs, zarr_url, version_prefix, parquet_data)
        for array_name, array_data in zarr_arrays.items():
            subfolder = "results" if "/" in array_name else "masks"
            save_zarr_array(fs, zarr_url, version_prefix, array_name, array_data, subfolder=subfolder)

        update_global_index(fs, zarr_url, version_prefix, user_id, metadata)

        return jsonify({
            'status': 'Version created',
            'zarr_url': zarr_url,
            'version_prefix': version_prefix,
            'user_id': user_id,
            'metadata_keys': list(metadata.keys()),
            'version_path': f"{zarr_url}/versions/{version_prefix}"
        })

    except Exception as e:
        app.logger.error(f"Error in add_version: {str(e)}")
        return jsonify({'error': str(e)}), 500

@app.route('/list_versions', methods=['GET'])
def list_versions():
    try:
        zarr_url = request.args.get('zarr_url')
        credential_id = request.args.get('credential_id', 'default')
        
        if not zarr_url:
            return jsonify({'error': 'Missing zarr_url parameter'}), 400

        fs = get_filesystem(credential_id)
        versions_dir = f"{zarr_url}/versions"
        
        if not fs.exists(versions_dir):
            app.logger.debug(f"No versions directory found: {versions_dir}")
            return jsonify({'versions': []})

        versions = []
        version_folders = [f['name'].split('/')[-1] for f in fs.ls(versions_dir, detail=True) if f['type'] == 'directory']
        app.logger.debug(f"Found version folders: {version_folders}")
        
        for version_prefix in version_folders:
            zattrs_path = f"{versions_dir}/{version_prefix}/.zattrs"
            try:
                if fs.exists(zattrs_path):
                    with fs.open(zattrs_path, 'r') as f:
                        metadata = json.load(f)
                    versions.append({
                        'version_prefix': version_prefix,
                        'metadata': metadata,
                        'path': f"{versions_dir}/{version_prefix}"
                    })
                else:
                    app.logger.warning(f"No .zattrs found for {version_prefix}")
            except Exception as e:
                app.logger.error(f"Error reading {zattrs_path}: {str(e)}")
        
        return jsonify({'versions': versions})
        
    except Exception as e:
        app.logger.error(f"Error in list_versions: {str(e)}")
        return jsonify({'error': str(e)}), 500

@app.route('/check_version_status', methods=['GET'])
def check_version_status():
    try:
        zarr_url = request.args.get('zarr_url')
        version_prefix = request.args.get('version_prefix')
        credential_id = request.args.get('credential_id', 'default')

        if not zarr_url or not version_prefix:
            return jsonify({'error': 'Missing zarr_url or version_prefix'}), 400

        fs = get_filesystem(credential_id)
        version_path = f"{zarr_url}/versions/{version_prefix}"
        zattrs_path = f"{version_path}/.zattrs"
        
        # 1. Check if the folder exists
        if not fs.exists(version_path):
            return jsonify({'status': 'missing', 'version_prefix': version_prefix}), 200

        # 2. Check if the metadata file exists and is readable
        metadata = None
        if fs.exists(zattrs_path):
            with fs.open(zattrs_path, 'r') as f:
                metadata = json.load(f)
            
            return jsonify({
                'status': 'verified',
                'version_prefix': version_prefix,
                'metadata': metadata,
                'message': 'Version folder and metadata confirmed.'
            }), 200
        else:
            return jsonify({
                'status': 'folder_exists_no_metadata',
                'version_prefix': version_prefix,
                'message': 'Version folder exists, but .zattrs metadata file is missing or unreadable.'
            }), 200

    except Exception as e:
        app.logger.error(f"Error checking version status: {str(e)}")
        return jsonify({'error': f"Internal check error: {str(e)}"}), 500

@app.route('/update_results', methods=['POST'])
def update_results():
    try:
        data = request.get_json()
        zarr_url = data.get('zarr_url')
        user_id = data.get('user_id')
        version_prefix = data.get('version_prefix')
        parquet_data = data.get('parquet_data', {})
        credential_id = data.get('credential_id', 'default')
        
        if not zarr_url or not version_prefix:
            return jsonify({'error': 'Missing zarr_url or version_prefix'}), 400
        
        fs = get_filesystem(credential_id)
        version_path = f"{zarr_url}/versions/{version_prefix}"
        if not fs.exists(version_path):
            return jsonify({'error': f"Version {version_prefix} does not exist"}), 404
        
        # Save Parquet files
        for key in parquet_data.keys():
            parquet_path = f"{version_path}/results/{key}.parquet"
            if fs.exists(parquet_path):
                fs.rm(parquet_path)
            df = pd.DataFrame(parquet_data[key])
            with fs.open(parquet_path, 'wb') as f:
                df.to_parquet(f)
        
        # Update metadata
        zattrs_path = f"{version_path}/.zattrs"
        metadata = data.get('metadata', {})
        with fs.open(zattrs_path, 'w') as f:
            json.dump(metadata, f)
        
        return jsonify({
            'status': 'Results updated',
            'version_prefix': version_prefix,
            'zarr_url': zarr_url,
            'user_id': user_id,
            'updated_parquet_keys': list(parquet_data.keys())
        })
    except Exception as e:
        app.logger.error(f"Error in update_results: {str(e)}")
        return jsonify({'error': str(e)}), 500

@app.route('/list_results', methods=['GET'])
def list_results():
    try:
        zarr_url = request.args.get('zarr_url')
        credential_id = request.args.get('credential_id', 'default')
        
        if not zarr_url:
            return jsonify({'error': 'Missing zarr_url parameter'}), 400

        fs = get_filesystem(credential_id)
        versions_dir = f"{zarr_url}/versions"
        
        if not fs.exists(versions_dir):
            app.logger.debug(f"No versions directory found: {versions_dir}")
            return jsonify({'results': []})

        results = []
        version_folders = [f['name'].split('/')[-1] for f in fs.ls(versions_dir, detail=True) if f['type'] == 'directory']
        app.logger.debug(f"Found version folders: {version_folders}")
        
        for version_prefix in version_folders:
            results_path = f"{versions_dir}/{version_prefix}/results"
            zattrs_path = f"{versions_dir}/{version_prefix}/.zattrs"
            try:
                metadata = {}
                if fs.exists(zattrs_path):
                    with fs.open(zattrs_path, 'r') as f:
                        metadata = json.load(f)
                
                parquet_files = []
                if fs.exists(results_path):
                    parquet_files = [
                        f['name'].split('/')[-1]
                        for f in fs.ls(results_path, detail=True)
                        if f['type'] == 'file' and f['name'].endswith('.parquet')
                    ]
                
                if parquet_files:
                    results.append({
                        'version_prefix': version_prefix,
                        'metadata': metadata,
                        'path': f"{versions_dir}/{version_prefix}",
                        'parquet_files': parquet_files
                    })
            except Exception as e:
                app.logger.error(f"Error processing version {version_prefix}: {str(e)}")
                continue
        
        return jsonify({'results': results})
        
    except Exception as e:
        app.logger.error(f"Error in list_results: {str(e)}")
        return jsonify({'error': str(e)}), 500

@app.route('/upload_parquet', methods=['POST'])
def upload_parquet():
    try:
        # Get form data
        zarr_url = request.form.get('zarr_url')
        user_id = request.form.get('user_id')
        version_prefix = request.form.get('version_prefix')
        parquet_key = request.form.get('parquet_key')
        credential_id = request.form.get('credential_id', 'default')
        metadata_json = request.form.get('metadata')  # 🆕 GET METADATA
        
        if not all([zarr_url, version_prefix, parquet_key]):
            return jsonify({'error': 'Missing required fields'}), 400
        
        # 🆕 Parse metadata
        metadata = {}
        if metadata_json:
            try:
                metadata = json.loads(metadata_json)
                app.logger.debug(f"Received metadata: {metadata}")
            except json.JSONDecodeError as e:
                app.logger.warning(f"Failed to parse metadata: {str(e)}")
        
        # Get uploaded file
        if 'file' not in request.files:
            return jsonify({'error': 'No file uploaded'}), 400
        
        file = request.files['file']
        if file.filename == '':
            return jsonify({'error': 'Empty filename'}), 400
        
        # Validate version exists
        fs = get_filesystem(credential_id)
        version_path = f"{zarr_url}/versions/{version_prefix}"
        if not fs.exists(version_path):
            return jsonify({'error': f"Version {version_prefix} does not exist"}), 404
        
        # Save to temporary file
        with tempfile.NamedTemporaryFile(delete=False, suffix='.parquet') as tmp:
            file.save(tmp.name)
            tmp.flush()
            
            # Get local file size for verification
            local_size = os.path.getsize(tmp.name)
            app.logger.debug(f"Local temp file size: {local_size} bytes")
            
            # Upload directly to S3
            results_dir = f"{version_path}/results"
            if not fs.exists(results_dir):
                fs.mkdir(results_dir)
            
            parquet_path = f"{results_dir}/{parquet_key}.parquet"
            
            # Copy file directly to S3
            with open(tmp.name, 'rb') as local_f:
                with fs.open(parquet_path, 'wb') as remote_f:
                    shutil.copyfileobj(local_f, remote_f)
            
            # Verify uploaded size
            remote_info = fs.info(parquet_path)
            remote_size = remote_info.get('size', 0)
            app.logger.debug(f"Remote S3 file size: {remote_size} bytes")
            
            # 🆕 SAVE METADATA to version's .zattrs
            if metadata:
                zattrs_path = f"{version_path}/.zattrs"
                
                # Read existing metadata
                existing_metadata = {}
                if fs.exists(zattrs_path):
                    with fs.open(zattrs_path, 'r') as f:
                        existing_metadata = json.load(f)
                
                # Update with new result metadata
                if 'dge_results' not in existing_metadata:
                    existing_metadata['dge_results'] = {}
                
                existing_metadata['dge_results'][parquet_key] = metadata
                
                # Write back
                with fs.open(zattrs_path, 'w') as f:
                    json.dump(existing_metadata, f, indent=2)
                
                app.logger.debug(f"Updated .zattrs with metadata for {parquet_key}")
            
            # Clean up temp file
            os.unlink(tmp.name)
            
            # Size verification
            if abs(remote_size - local_size) > 1024:  # Allow 1KB difference
                app.logger.warning(f"Size mismatch: local={local_size}, remote={remote_size}")
            
            return jsonify({
                'status': 'Results uploaded',
                'version_prefix': version_prefix,
                'zarr_url': zarr_url,
                'user_id': user_id,
                'parquet_key': parquet_key,
                'uploaded_size_kb': round(remote_size / 1024, 2),
                'local_size_kb': round(local_size / 1024, 2),
                'metadata_saved': bool(metadata)  # 🆕 CONFIRM METADATA SAVED
            })
    
    except Exception as e:
        app.logger.error(f"Error in upload_parquet: {str(e)}")
        return jsonify({'error': str(e)}), 500


@app.route('/get_result_metadata', methods=['GET'])
def get_result_metadata():
    try:
        zarr_url = request.args.get('zarr_url')
        version_prefix = request.args.get('version_prefix')
        parquet_key = request.args.get('parquet_key')
        credential_id = request.args.get('credential_id', 'default')
        
        if not all([zarr_url, version_prefix, parquet_key]):
            return jsonify({'error': 'Missing required parameters'}), 400
        
        fs = get_filesystem(credential_id)
        
        # Read version metadata
        zattrs_path = f"{zarr_url}/versions/{version_prefix}/.zattrs"
        if not fs.exists(zattrs_path):
            return jsonify({'error': 'Metadata not found'}), 404
        
        with fs.open(zattrs_path, 'r') as f:
            version_metadata = json.load(f)
        
        # 🆕 GET RESULT-SPECIFIC METADATA
        result_metadata = version_metadata.get('dge_results', {}).get(parquet_key, {})
        
        # Get parquet file info
        parquet_path = f"{zarr_url}/versions/{version_prefix}/results/{parquet_key}.parquet"
        if fs.exists(parquet_path):
            file_info = fs.info(parquet_path)
            result_metadata['file_size_kb'] = round(file_info.get('size', 0) / 1024, 2)
            
            # Optional: Read parquet to get gene count if not in metadata
            if 'n_genes' not in result_metadata:
                try:
                    with fs.open(parquet_path, 'rb') as f:
                        df = pd.read_parquet(f)
                        result_metadata['n_genes'] = len(df)
                        if 'pvals_adj' in df.columns:
                            result_metadata['n_significant'] = int((df['pvals_adj'] < 0.05).sum())
                except Exception as e:
                    app.logger.warning(f"Could not read parquet: {str(e)}")
        
        return jsonify(result_metadata)
        
    except Exception as e:
        app.logger.error(f"Error in get_result_metadata: {str(e)}")
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8080, debug=True)