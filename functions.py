import sqlite3
import pandas as pd
import json

def getPandasDFfromDB(pathToDBfile):
	#My friend Murat wrote these two lines, I have no idea what the squlite3 routine is doing
	conn = sqlite3.connect(pathToDBfile)
	database_df = pd.read_sql("SELECT * FROM calculations", conn)
	return database_df

def getXYZstructureList(pathToDBfile):
	df = getPandasDFfromDB(pathToDBfile)
	optimized_xyz_list = []
	for index, row in df.iterrows():
		try:
			blob_data_dict = json.loads(row['blob_data'])
			optimized_xyz_list.append(blob_data_dict['opt_xyz'])
		except (json.JSONDecodeError, KeyError) as  e:
			print(f"Error processing row {index}: {e}")
			optimized_xyz_list.append(None)
	return optimized_xyz_list

