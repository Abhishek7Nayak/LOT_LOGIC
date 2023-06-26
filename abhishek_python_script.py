import pandas as pd
from datetime import datetime, timedelta

# Read the CSV file into a pandas DataFrame
df = pd.read_csv("/home/toy_dataset.csv")

# Convert date columns to datetime format
df["drug_start"] = pd.to_datetime(df["drug_start"])
df["drug_end"] = pd.to_datetime(df["drug_end"])

# Get the minimum drug start date for each patient
df_patient_min_start = df.groupby('patid')['drug_start'].min().reset_index()

# Create an empty DataFrame to store the transformed data
df_transformed = pd.DataFrame(columns=['patid', 'treatment', 'drug_start', 'drug_end', 'lot_new'])

# Process each patient
for i in range(len(df_patient_min_start)):
    patient_id = df_patient_min_start.loc[i, 'patid']
    min_start_date = df_patient_min_start.loc[i, 'drug_start']

    # Filter treatments within 30 days from the first treatment start date
    df_30_days = df[(df['patid'] == patient_id) & (df['drug_start'] <= (min_start_date + timedelta(days=30)))].reset_index()

    # Assign lot 1 for the treatments within 30 days
    df_30_days['lot_new'] = 1

    # Find unique treatments within the 30-day window
    unique_treatments = set()
    for _, row in df_30_days.iterrows():
        treatments = row['treatment'].split(' + ')
        unique_treatments.update(treatments)

    # Create a string representation of the unique treatments
    unique_treatments_string = ' + '.join(unique_treatments)

    # Create a new DataFrame for the line 1 data
    df_line1 = pd.DataFrame({
        'patid': [patient_id],
        'treatment': unique_treatments_string,
        'drug_start': df_30_days['drug_start'].min(),
        'drug_end': df_30_days['drug_end'].max(),
        'lot_new': df_30_days['lot_new'].unique()
    })

    # Update the transformed DataFrame with line 1 data
    df_transformed = pd.concat([df_transformed, df_line1])

    # Filter treatments after 30 days
    df_after_30_days = df[(df['patid'] == patient_id) & (df['drug_start'] > df_line1['drug_end'].iloc[0])].reset_index()

    # Initialize lot count
    lot_count = 1

    # Process each treatment after 30 days
    for _, row in df_after_30_days.iterrows():
        treatments = row['treatment'].split(" + ")

        # Check if any treatment is not present in the previous treatments
        if not any(treatment in unique_treatments for treatment in treatments):
            lot_count += 1

        # Update the lot number for the treatment
        df_after_30_days.at[_, 'lot_new'] = lot_count

        # Update the unique treatments set
        unique_treatments.update(treatments)

    # Update the transformed DataFrame with the treatments after 30 days
    df_transformed = pd.concat([df_transformed, df_after_30_days])

# Fill null values in the lot_new column with 1
df_transformed['lot_new'].fillna(1, inplace=True)

# Select required columns and display the transformed DataFrame
df_transformed = df_transformed[['patid', 'treatment', 'drug_start', 'drug_end', 'lot_new']]
#print(df_transformed.head())

#
# Save the transformed DataFrame to a CSV file
df_transformed.to_csv('lotfile_all.csv', index=False)
