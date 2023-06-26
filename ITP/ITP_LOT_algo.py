import pandas as pd

# Read the input data from a Parquet file
df = pd.read_parquet("/home/test.parquet")

# Define variable names
cor = 'corticosteroids'
ig = 'immunoglobulin'
treatment_col = 'treatment'
pat_id = 'hvid'
lot_col = 'new_lot'
TPO = ['avatrombopag', 'eltrombopag', 'romiplostim']
rituximab = 'rituximab'

# Sort the DataFrame by patient ID and start date
df = df.sort_values(by=[pat_id, 'start_date'])
final = pd.DataFrame()

# Process each unique patient ID
for patient_id in df[pat_id].unique():
    # Extract the data for the current patient
    df_patient = df[df[pat_id] == patient_id]
    df_patient = df_patient.reset_index()
    df_patient[lot_col] = 9999999
    RT = 888888

    initial_lot = 1
    list_store_treatment = []

    # Iterate over each row for the current patient
    for index, row in df_patient.iterrows():
        current_treatment = row[treatment_col]

        # Check if the current treatment contains corticosteroids, immunoglobulin, or their combinations
        if any(treatment in current_treatment for treatment in [cor, ig, 'immunoglobulin + corticosteroids',
                                                               'corticosteroids + immunoglobulin',
                                                               'corticosteroids + corticosteroids',
                                                               'corticosteroids + corticosteroids + immunoglobulin',
                                                               'corticosteroids + immunoglobulin + corticosteroids',
                                                               'antid', 'antid + corticosteroids',
                                                               'corticosteroids + antid']):
            # Assign initial lot number for the first occurrence of corticosteroids or immunoglobulin
            if initial_lot == 1:
                df_patient.loc[index, lot_col] = initial_lot
                list_store_treatment.append(current_treatment)
                df_patient.loc[index, "new_treat"] = current_treatment
            # Assign RT (Rescue Therapy) lot number for subsequent occurrences
            elif df_patient.loc[index - 1, lot_col] > 1:
                df_patient.loc[index, lot_col] = RT
                list_store_treatment.append(current_treatment)
                df_patient.loc[index, "new_treat"] = current_treatment
            # Assign initial lot number for subsequent occurrences of corticosteroids if not combined with immunoglobulin
            elif current_treatment == cor and not any(treatment in list_store_treatment for treatment in [cor, ig]):
                df_patient.loc[index, lot_col] = initial_lot
                list_store_treatment.append(current_treatment)
                df_patient.loc[index, "new_treat"] = current_treatment

        # Check if the current treatment contains eltrombopag
        elif TPO[1] in current_treatment.split(" + "):
            # Assign new lot number for the first occurrence of eltrombopag
            if row[lot_col] == 1:
                initial_lot += 1
                df_patient.loc[index, lot_col] = initial_lot
                list_store_treatment.append(TPO[1])
                df_patient.loc[index, "new_treat"] = TPO[1]
            # Assign RT (Rescue Therapy) lot number for subsequent occurrences
            elif row[lot_col] > 1 and df_patient.loc[index - 1, "new_treat"] == TPO[1]:
                df_patient.loc[index, lot_col] = RT
                list_store_treatment.append(TPO[1])
                df_patient.loc[index, "new_treat"] = TPO[1]

        # Check if the current treatment contains avatrombopag
        elif TPO[0] in current_treatment.split(" + "):
            # Assign new lot number for the first occurrence of avatrombopag
            if row[lot_col] == 1:
                initial_lot += 1
                df_patient.loc[index, lot_col] = initial_lot
                list_store_treatment.append(TPO[0])
                df_patient.loc[index, "new_treat"] = TPO[0]
            # Assign RT (Rescue Therapy) lot number for subsequent occurrences
            elif row[lot_col] > 1 and df_patient.loc[index - 1, "new_treat"] == TPO[0]:
                df_patient.loc[index, lot_col] = RT
                list_store_treatment.append(TPO[0])
                df_patient.loc[index, "new_treat"] = TPO[0]

        # Check if the current treatment contains romiplostim
        elif TPO[2] in current_treatment.split(" + "):
            # Assign new lot number for the first occurrence of romiplostim
            if row[lot_col] == 1:
                initial_lot += 1
                df_patient.loc[index, lot_col] = initial_lot
                list_store_treatment.append(TPO[2])
                df_patient.loc[index, "new_treat"] = TPO[2]
            # Assign RT (Rescue Therapy) lot number for subsequent occurrences
            elif row[lot_col] > 1 and df_patient.loc[index - 1, "new_treat"] == TPO[2]:
                df_patient.loc[index, lot_col] = RT
                list_store_treatment.append(TPO[2])
                df_patient.loc[index, "new_treat"] = TPO[2]

        elif rituximab in current_treatment.split(" + "):
            # Assign new lot number for the first occurrence of rituximab
            if row[lot_col] == 1:
                initial_lot += 1
                df_patient.loc[index, lot_col] = initial_lot
                list_store_treatment.append(rituximab)
                df_patient.loc[index, "new_treat"] = rituximab

    # Concatenate the processed data for the current patient with the final DataFrame
    final = pd.concat([final, df_patient])

# Remove the 'index' column from the final DataFrame
final = final.drop(columns=['index'])

# Perform post-processing
t = final
final_postproc = pd.DataFrame()

# Process each unique patient ID
for i in t[pat_id].unique():
    t1 = t[t[pat_id] == i]
    t1 = t1.reset_index()
    list_store_treatment = []
    # Iterate over each row for the current patient
    for j in range(t1.shape[0]):
        if t1.loc[j, lot_col] != 888888:
            t1.loc[j, 'test_lot'] = t1.loc[j, lot_col]
            list_store_treatment.append(t1.loc[j, lot_col])
        else:
            t1.loc[j, 'test_lot'] = str(list_store_treatment[-1]) + "_" + "RT"
    final_postproc = t1.append(final_postproc, ignore_index=True)

# Select the desired columns and rename them
t2 = final_postproc[['hvid', 'treatment', 'start_date', 'drug_end', 'new_lot', 'test_lot', 'new_treat']]
t2.rename(columns={'new_lot': 'lot_int', 'test_lot': 'lot'}, inplace=True)

# Save the final DataFrame to a CSV file
t2.to_csv('lotfile_all.csv', index=False)
