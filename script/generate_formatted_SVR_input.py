import numpy as np

MAX_RADIUS = 55

def get_top_n_features(threshold, feature_ranks, X):
    '''
    This method gets the top 'threshold' ranked features according to their Pearson Pearson

    Parameters:
    -----------
    threshold: int
        This is the top n to retrieve based on correlation

    feature_ranks: list[string]
        this is the ranks(ordered best-> worst) of the features from the flattened matrix

    X: list[np.ndarray]
        This is the input data, a list of flattened numpy arrays

    Return:
    -----------
    list[np.ndarray]
        This method returns a list of 1D numpy arrays that have the correct ranking top n features

    '''
    feature_training_set = []
    feature_testing_set = []

    for training_example in X:
        augmented_example = []
        for feature_number in feature_ranks[:threshold]:
            augmented_example.append(training_example[int(feature_number)])

        feature_training_set.append(np.asarray(augmented_example))


    return feature_training_set

def get_feature_ranks():
    '''
    This method returns the ordered feature ranks, precomputed and saved

    Parameters:
    -------------
    pathToFeatureScores: string
        This is a hard-coded path, points to the Pearson_Correlation_Individula_Features.txt

    Returns:
    -----------
    list: [string]
        A list of strings of the feature indexes ranked from best to worst by pearson correlation
    '''
    #load and parse the feature ranks
    pathToFeatureScores = './script/Pearson_Correlation_Individula_Features.txt'

    raw_data = open(pathToFeatureScores).read()
    feature_ranks = []
    for line in raw_data.split('\n'):
        line_data = line.split('\t')
        feature_ranks.append(line_data[0])

    return feature_ranks

def flatten(data):
    '''
    This method flattens each input data into a vector

    Parameters:
    ----------
    data: list[np.ndarray((47,51))]
        This is a list of the input examples, 47 rows, 51 columns

    Return:
    ------
    This method does not return anything, it modifies the data parameter by reference (save on RAM)
    '''
    for i in range(len(data)):
        x_flatten = data[i].flatten()
        data[i] = x_flatten

def parse_server_data(server_data, top_n):
    '''
    This method is responsible for creating the final training and label data.
    Shape for the input data: (47,51)

    Parameters:
    -----------
    server_data: dictionary
        This is one file created by the 'generate_casp_fragment_structures.py' script. It is all the relevant data for a server prediction for a
        target

    Return:
    list[np.ndarray(47,51)], list[float]
        this method returns two items. The first is a list of 51x21 matrix reperesenting the input data for each residue in the target
        and the second is the localQA score for the corresponding input data
    '''
    feature_ranks = get_feature_ranks()

    server_X,server_y = [],[]

    for index, data_dictionary in server_data.items():
        X = []
        y = data_dictionary['local_qa']

        density_change_data = np.transpose(data_dictionary['aa_density_change'])
        for aa_row in density_change_data:
            X.append(aa_row[:MAX_RADIUS])

        hydro_change = data_dictionary['hydro_change']
        X.append(hydro_change[:MAX_RADIUS])

        mass_change = data_dictionary['mass_change']
        X.append(mass_change[:MAX_RADIUS])

        sol_change = data_dictionary['sol_change']
        X.append(sol_change[:MAX_RADIUS])

        iso_change = data_dictionary['iso_change']
        X.append(iso_change[:MAX_RADIUS])

        ave_dist = data_dictionary['average_distance']
        X.append(ave_dist[:MAX_RADIUS])

        std_dist = data_dictionary['std_dev_distance']
        X.append(std_dist[:MAX_RADIUS])

        percent_contact = data_dictionary['percent_contact']
        X.append(percent_contact[:MAX_RADIUS])

        contact_fequency_map = np.transpose(data_dictionary['structure_contact_matrix'])
        for aa_row in contact_fequency_map:
            X.append(aa_row[:MAX_RADIUS])

        #misc rows
        #total entries: 4 + 1 + 1 + 1 + 1 + 1 + 1 + 3 + 20 = 33
        rf_pred = data_dictionary['rf_predictions']
        aa_mass, aa_hydro, aa_sol, aa_iso = data_dictionary['aa_mass'], data_dictionary['aa_hydro'], data_dictionary['aa_sol'], data_dictionary['aa_iso']
        psi, phi = data_dictionary['psiphi'][0], data_dictionary['psiphi'][1]
        ss = data_dictionary['ss_encoded']
        aa_one_hot = data_dictionary['aa_encoded']

        row_1 = []
        row_1.extend(rf_pred)
        row_1.extend([aa_mass, aa_hydro, aa_sol,])
        row_1.extend(aa_iso)
        row_1.extend([psi, phi])
        row_1.extend(ss)
        row_1.extend([0]* (MAX_RADIUS - len(row_1)))

        row_2 = []
        row_2.extend(aa_one_hot)
        row_2.extend([0]* (MAX_RADIUS - len(row_2)))

        server_X.append(np.asarray(X))
        server_y.append(np.asarray(y))

    flatten(server_X)
    svr_input = get_top_n_features(top_n,feature_ranks,  server_X)

    return np.asarray(svr_input), np.array(server_y)
