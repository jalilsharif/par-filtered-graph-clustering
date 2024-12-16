#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// Function to generate a symmetrical matrix and a dissimilarity matrix
void generateMatrices(const string &inputFile, const string &symMatrixFile, const string &dissMatrixFile, const string &metadataFile) {
    ifstream input(inputFile);
    if (!input.is_open()) {
        cerr << "Error: Cannot open input file " << inputFile << endl;
        exit(1);
    }

    unordered_map<string, int> node_to_index;
    vector<tuple<int, int, double>> edges;
    string row, col;
    double rho;
    int index = 0;

    // Read adjacency list and map node labels to indices
    string header;
    getline(input, header); // Skip header row
    while (input >> row >> col >> rho) {
        if (node_to_index.find(row) == node_to_index.end()) {
            node_to_index[row] = index++;
        }
        if (node_to_index.find(col) == node_to_index.end()) {
            node_to_index[col] = index++;
        }
        edges.emplace_back(node_to_index[row], node_to_index[col], rho);
    }
    input.close();

    size_t n = node_to_index.size();

    // Initialize symmetrical matrix and dissimilarity matrix
    vector<vector<double>> symMatrix(n, vector<double>(n, 0.0));
    vector<vector<double>> dissMatrix(n, vector<double>(n, 0.0));

    for (const auto &[i, j, weight] : edges) {
        symMatrix[i][j] = weight;
        symMatrix[j][i] = weight;

        double diss_value = sqrt(2 * (1.0 - weight));
        dissMatrix[i][j] = diss_value;
        dissMatrix[j][i] = diss_value;
    }
// Set diagonal elements to 1 for symmetrical matrix
    for (size_t i = 0; i < n; ++i) {
        symMatrix[i][i] = 1.0;
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

    metaFile << "Unique_IDS " << n << endl;
    for (const auto &[label, idx] : node_to_index) {
        metaFile << (idx+1) << " " << label << endl;
    }
    metaFile.close();

    cout << "Matrix generation completed successfully." << endl;
    cout << "Symmetrical matrix saved to: " << symMatrixFile << endl;
    cout << "Dissimilarity matrix saved to: " << dissMatrixFile << endl;
    cout << "Metadata saved to: " << metadataFile << endl;
}

int main(int argc, char *argv[]) {
    if (argc < 5) {
        cerr << "Usage: " << argv[0] << " <input_file> <sym_matrix_file> <diss_matrix_file> <metadata_file>" << endl;
        return 1;
    }

    string inputFile = argv[1];
    string symMatrixFile = argv[2];
    string dissMatrixFile = argv[3];
    string metadataFile = argv[4];

    generateMatrices(inputFile, symMatrixFile, dissMatrixFile, metadataFile);

    return 0;
}
