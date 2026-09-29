import os
import pandas as pd

global_dataframes = {}
relevant_csv_files = {
    "cell_outline",
    "Golgi_convex_hull",
    "nuclei",
    "golgi_stack",
    "golgi_disconnected",
}
channels = ["KO", "HAP40KO", "Q23", "Q145"]
batch_id_counter = {"KO": 0, "HAP40KO": 0, "Q23": 0, "Q145": 0}


def update_global_dataframes(df, key, batch_id, channel_id):
    batch = f"{channel_id}_{batch_id}"
    df["batch_id"] = batch
    df["channel_id"] = channel_id

    if key not in global_dataframes:
        global_dataframes[key] = df
    else:
        global_dataframes[key] = pd.concat(
            [global_dataframes[key], df], ignore_index=True
        )
    print(
        f"Updated global dataframes for {key} with {len(df)} additional columns from experiment #{batch_id}, current length: {len(global_dataframes[key])}"
    )
    return global_dataframes


def read_csv_files(directory):
    """
    Read all csv files in a directory and return a dictionary of dataframes
    :param directory: directory path
    :return: dataframes dictionary including all csv files in the directory
    """
    global batch_id_counter
    dataframes = {}
    channel_id = "unknown"
    for channel in channels:
        if channel.lower() in directory.lower():
            channel_id = channel

    batch_id = batch_id_counter[channel_id]

    csv_files = [file for file in os.listdir(directory) if file.endswith(".csv")]

    for file in csv_files:
        file_path = os.path.join(directory, file)
        dataframe_name = os.path.splitext(file)[0]
        if dataframe_name not in relevant_csv_files:
            print(f"Skipping {dataframe_name}")
            continue
        dataframes[dataframe_name] = pd.read_csv(file_path)
        update_global_dataframes(
            dataframes[dataframe_name], dataframe_name, batch_id, channel_id
        )
    batch_id_counter[channel_id] += 1

    return dataframes


def calculate(dfs, directory_path, object_identifiers, channel_id=None):
    """
    Calculate the analysis for a given directory
    :param dfs: all dataframes in the directory
    :param directory_path: path to write the results to
    :return: dataframe including the analysis
    """
    print(f"Starting calculation for {directory_path}")

    convex_hull = dfs["Golgi_convex_hull"][
        object_identifiers + ["AreaShape_Area"]
    ].rename(columns={"AreaShape_Area": "convex_hull_area"})

    nuclei = dfs["nuclei"][object_identifiers + ["AreaShape_Area"]].rename(
        columns={"AreaShape_Area": "nucleus_area"}
    )

    golgi_stack = dfs["golgi_disconnected"][
        object_identifiers + ["AreaShape_Area"]
    ].rename(columns={"AreaShape_Area": "golgi_stack_area"})

    object_identifiers2 = object_identifiers.copy()
    object_identifiers2.remove("ObjectNumber")
    golgi_fragments = (
        dfs["golgi_stack"][object_identifiers + ["Parent_cell_outline"]]
        .groupby(object_identifiers2 + ["Parent_cell_outline"])
        .count()
        .reset_index()
        .rename(
            columns={
                "ObjectNumber": "fragments_count",
                "Parent_cell_outline": "ObjectNumber",
            }
        )
    )

    result = (
        dfs["cell_outline"]
        .merge(convex_hull, on=object_identifiers, how="left")
        .merge(nuclei, on=object_identifiers, how="left")
        .merge(golgi_stack, on=object_identifiers, how="left")
        .merge(golgi_fragments, on=object_identifiers, how="left")
    )

    result["golgi_area"] = result["convex_hull_area"] / result["nucleus_area"]
    result["golgi_compactness"] = (
        result["golgi_stack_area"] / result["convex_hull_area"]
    )

    result.rename(
        columns={
            "ObjectNumber": "cell_number",
            "ImageNumber": "image_number",
        }
    )

    if channel_id:
        print(f"Filtering for channel: {channel_id}")
        result = result[result["channel_id"] == channel_id]

    col_to_calc = ["golgi_area", "golgi_compactness", "fragments_count"]

    for selected_column in col_to_calc:
        print(f"Calculating z-score for {selected_column}")
        mean_value = result[selected_column].mean()
        std_deviation = result[selected_column].std()
        print(f"Mean: {mean_value}")
        print(f"Standard deviation: {std_deviation}")
        print("=============")

        result[f"{selected_column}_zscore"] = (
            result[selected_column] - mean_value
        ) / std_deviation

    range_threshold = 2

    col_to_calc_2 = [f"{selected_column}_zscore" for selected_column in col_to_calc]

    result["in_range"] = result[col_to_calc_2].apply(
        lambda row: "no" if any(abs(val) > range_threshold for val in row) else "yes",
        axis=1,
    )

    result_in_range = result[result["in_range"] == "yes"]
    out_of_range_df = result[result["in_range"] == "no"]

    if channel_id:
        path_in_range = os.path.join(
            directory_path, f"{channel_id}_analysis_in_range_.csv"
        )
        path_out_of_range = os.path.join(
            directory_path, f"{channel_id}_analysis_out_of_range_.csv"
        )
        path_with_zscore = os.path.join(
            directory_path, f"{channel_id}_analysis_with_zscore_.csv"
        )
    else:
        path_in_range = os.path.join(directory_path, "analysis_in_range.csv")
        path_out_of_range = os.path.join(directory_path, "analysis_out_of_range.csv")
        path_with_zscore = os.path.join(directory_path, "analysis_with_zscore.csv")

    result_in_range.to_csv(path_in_range, index=False)
    out_of_range_df.to_csv(path_out_of_range, index=False)
    result.to_csv(path_with_zscore, index=False)

    print("====================")
    print("====================")

    return result


def dfs_search_directory(directory_path):
    """
    Search for all subdirectories in a given directory until it finds a csv file
    :param directory_path: starting directory path
    :return: yields a directory path
    """
    for root, dirs, files in os.walk(directory_path):
        for file in files:
            if file.endswith(".csv") and (
                os.path.splitext(file)[0] in relevant_csv_files
            ):
                yield root
                break


def main():
    print(f"\nStarting search and analysis for sub-directories\n")
    for directory_path in dfs_search_directory("./"):
        print(directory_path)
        dfs = read_csv_files(directory_path)
        calculate(
            dfs, directory_path, object_identifiers=["ImageNumber", "ObjectNumber"]
        )
    print(f"\nStarting analysis for all files combined\n")

    for channel in channels:
        calculate(
            global_dataframes,
            "./",
            object_identifiers=["ImageNumber", "ObjectNumber", "batch_id"],
            channel_id=channel,
        )


if __name__ == "__main__":
    main()
