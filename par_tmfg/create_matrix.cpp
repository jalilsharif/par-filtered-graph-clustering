#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>

using namespace std;

// Helper function to append "_dissimilarity" to the file name
string generateDissimilarityFileName(const string &symMatrixFile) {
    size_t dotPos = symMatrixFile.find_last_of(".");
    if (dotPos != string::npos) {
        return symMatrixFile.substr(0, dotPos) + "_dissimilarity" + symMatrixFile.substr(dotPos);
    } else {
        return symMatrixFile + "_dissimilarity";
    }
}

// Function to generate a symmetrical matrix and a dissimilarity matrix
void generateMatrices(const string &inputFile, const string &symMatrixFile, const string &dissMatrixFile, const string &metadataFile, int edgeLimit = -1) {
    ifstream input(inputFile);
    if (!input.is_open()) {
        cerr << "Error: Cannot open input file " << inputFile << endl;
        exit(1);
    }

    unordered_map<string, int> node_to_index;
    unordered_map<int, string> index_to_node; // Reverse mapping
    vector<tuple<int, int, double>> edges;
    string row, col;
    double rho;
    int index = 0;

    // Read adjacency list and map node labels to indices
    string header;
    getline(input, header); // Skip header row
    while (input >> row >> col >> rho) {
        if (node_to_index.find(row) == node_to_index.end()) {
            node_to_index[row] = index;
            index_to_node[index] = row; // Reverse mapping
            index++;
        }
        if (node_to_index.find(col) == node_to_index.end()) {
            node_to_index[col] = index;
            index_to_node[index] = col; // Reverse mapping
            index++;
        }
        edges.emplace_back(node_to_index[row], node_to_index[col], rho);
    }
    input.close();

    size_t n = node_to_index.size();

    // Initialize symmetrical matrix and dissimilarity matrix
    vector<vector<double>> symMatrix(n, vector<double>(n, 0.0));
    vector<vector<double>> dissMatrix(n, vector<double>(n, 0.0));

    // Process edges
    cout << "Edges with node labels and indices:" << endl;
    int printedEdges = 0;
    for (const auto &[i, j, weight] : edges) {
        if (edgeLimit > 0 && printedEdges >= edgeLimit) {
            break; // Stop printing if limit is reached
        }

        cout << "Edge between node " << i << " (" << index_to_node[i] << ") and node " << j << " (" << index_to_node[j] << ") with weight " << weight << endl;

        symMatrix[i][j] = weight;
        symMatrix[j][i] = weight;

        double diss_value = sqrt(2 * (1.0 - weight));
        dissMatrix[i][j] = diss_value;
        dissMatrix[j][i] = diss_value;

        printedEdges++;
    }



    // Write symmetrical matrix to binary file
    ofstream symFile(symMatrixFile, ios::out | ios::binary);
    if (!symFile.is_open()) {
        cerr << "Error: Cannot open symmetrical matrix output file " << symMatrixFile << endl;
        exit(1);
    }
    for (size_t i = 0; i < n; ++i) {
        symFile.write(reinterpret_cast<const char *>(symMatrix[i].data()), n * sizeof(double));
    }
    symFile.close();

    // Write dissimilarity matrix to binary file
    ofstream dissFile(dissMatrixFile, ios::out | ios::binary);
    if (!dissFile.is_open()) {
        cerr << "Error: Cannot open dissimilarity matrix output file " << dissMatrixFile << endl;
        exit(1);
    }
    for (size_t i = 0; i < n; ++i) {
        dissFile.write(reinterpret_cast<const char *>(dissMatrix[i].data()), n * sizeof(double));
    }
    dissFile.close();

    // Write metadata to file
    ofstream metaFile(metadataFile);
    if (!metaFile.is_open()) {
        cerr << "Error: Cannot open metadata file " << metadataFile << endl;
        exit(1);
    }

    metaFile << "Unique_IDs " << n << endl;
    for (const auto &[idx, label] : index_to_node) {
        metaFile << (idx+1) << " " << label << endl;
    }
    metaFile.close();

    cout << "Matrix generation completed successfully." << endl;
    cout << "Symmetrical matrix saved to: " << symMatrixFile << endl;
    cout << "Dissimilarity matrix saved to: " << dissMatrixFile << endl;
    cout << "Metadata saved to: " << metadataFile << endl;
}

int main(int argc, char *argv[]) {
    string inputFile;
    string symMatrixFile;
    string metadataFile;
    int edgeLimit = -1;

    // Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-i" && i + 1 < argc) {
            inputFile = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            symMatrixFile = argv[++i];
        } else if (arg == "-m" && i + 1 < argc) {
            metadataFile = argv[++i];
        } else if (arg == "-l" && i + 1 < argc) {
            edgeLimit = stoi(argv[++i]);
        } else {
            cerr << "Usage: " << argv[0]
                 << " -i <input_file> -o <sym_matrix_file> -m <metadata_file> [-l <edge_limit>]"
                 << endl;
            return 1;
        }
    }

    if (inputFile.empty() || symMatrixFile.empty() || metadataFile.empty()) {
        cerr << "Error: Missing required arguments." << endl;
        cerr << "Usage: " << argv[0]
             << " -i <input_file> -o <sym_matrix_file> -m <metadata_file> [-l <edge_limit>]"
             << endl;
        return 1;
    }

    // Generate dissimilarity matrix file name
    string dissMatrixFile = generateDissimilarityFileName(symMatrixFile);

    // Generate matrices
    generateMatrices(inputFile, symMatrixFile, dissMatrixFile, metadataFile, edgeLimit);

    return 0;
}
