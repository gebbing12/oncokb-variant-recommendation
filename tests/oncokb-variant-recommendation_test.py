import pytest
import pandas as pd
import boto3
import os
from botocore.exceptions import ClientError
import logging
import io


##test data
@pytest.fixture(params=[
    ("msk_impact", 
     "path_to_filter_file.txt", 
     False ) 
])

def test_data(request):
    '''
    param:
        test_name:test data sample name
        filter_path: Filter file path
        s3: whether to use s3
    '''
    return request.param

# Execlude data by SAMPLE ID
def exclude_samples(data,filter_path):
    excluded_sample_ids = pd.read_csv(filter_path,delimiter='\t',header=None)[0].tolist()
    filtered_data = data[~data['SAMPLE_ID'].isin(excluded_sample_ids)]
    return filtered_data
    

def dataframe_to_s3(dataframe, bucket, object_name):
    # Create a buffer
    csv_buffer = io.StringIO()
    dataframe.to_csv(csv_buffer, sep='\t', index=False)
    csv_buffer.seek(0)  

    # Upload to s3
    s3_client = boto3.client('s3')
    try:
        s3_client.put_object(Bucket=bucket, Key=object_name, Body=csv_buffer.getvalue())
        return True
    except Exception as e:
        logging.error(f"Failed to upload {object_name} to {bucket}. Error: {e}")
        return False



def analysis(test_name, filter_path=None,s3 = False):
    current_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(current_dir, '..','biodata/')
    #test_name="msk_impact_2017"
    file_path = path+test_name
    
    # get the data
    sample = pd.read_csv(file_path+"/data_clinical_sample.txt", usecols=['SAMPLE_ID', 'CANCER_TYPE','CANCER_TYPE_DETAILED'], delimiter='\t', skiprows=4)
    mutation = pd.read_csv(file_path+"/data_mutations.txt", usecols=['Hugo_Symbol', 'Tumor_Sample_Barcode','HGVSp_Short'], delimiter='\t',skiprows=1)
    data = pd.merge(sample, mutation, left_on='SAMPLE_ID', right_on='Tumor_Sample_Barcode').drop(columns=['Tumor_Sample_Barcode'])
    if filter_path:
        data = exclude_samples(data,filter_path)
        
    
    #all cancer frequency
    all_number = data[['Hugo_Symbol', 'HGVSp_Short']].value_counts().reset_index(name='Number')
    all_percentage = data[['Hugo_Symbol', 'HGVSp_Short']].value_counts(normalize=True).reset_index(name='Percentage')
    all_frequency = pd.merge(all_number, all_percentage, on=['Hugo_Symbol', 'HGVSp_Short'])
    
    #cancer type frequency
    type_number = data.groupby('CANCER_TYPE')[['Hugo_Symbol', 'HGVSp_Short']].value_counts().rename("Number").reset_index()
    type_percentage = data.groupby('CANCER_TYPE')[['Hugo_Symbol', 'HGVSp_Short']].value_counts(normalize=True).rename("Percentage")
    type_frequency = pd.merge(type_number, type_percentage, on=['CANCER_TYPE', 'Hugo_Symbol', 'HGVSp_Short'])
    
    #cancer type detailed frequency
    detailed_number = data.groupby(['CANCER_TYPE','CANCER_TYPE_DETAILED'])[['Hugo_Symbol', 'HGVSp_Short']].value_counts().rename("Number").reset_index()
    detailed_percentage = data.groupby(['CANCER_TYPE','CANCER_TYPE_DETAILED'])[['Hugo_Symbol', 'HGVSp_Short']].value_counts(normalize=True).rename("Percentage")
    detailed_frequency = pd.merge(detailed_number, detailed_percentage, on=['CANCER_TYPE', 'CANCER_TYPE_DETAILED', 'Hugo_Symbol', 'HGVSp_Short'])
    
    
    #save the data
    if s3:
        #upload to S3
        bucket = 'bucket_name'
        dataframe_to_s3(all_frequency, bucket, test_name + "_all_frequency.txt")
        dataframe_to_s3(type_frequency, bucket, test_name + "_type_frequency.txt")
        dataframe_to_s3(detailed_frequency, bucket, test_name + "_detailed_frequency.txt")
    else:
        # save to local
        save_path = os.path.join(current_dir, '..', 'test-data/') 
        detailed_frequency.to_csv(save_path + test_name + "_detailed_frequency.txt", sep='\t',index = False) 
        all_frequency.to_csv(save_path + test_name + "_all_frequency.txt", sep='\t',index = False) 
        type_frequency.to_csv(save_path + test_name + "_type_frequency.txt", sep='\t',index = False) 
    
    return True
        
def test_analysis(test_data):
    test_name, filter_path, s3 = test_data
    analysis(test_name, filter_path, s3)
    