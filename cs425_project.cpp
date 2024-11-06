#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>

using namespace std;

// 1. Apply Laplacian Operator to the "building.raw" and find all zero-crossing points.
// Function to apply the zero-crossing points using Laplacian operator
void applyLaplacian(const vector<unsigned char>& input, vector<unsigned char>& output, int width, int height, int threshold) {
    // Ensure the output vector has the correct size
    output.resize(input.size());

    // Loop through the image pixels, excluding the border pixels
    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            // Store pixel values in a 2x2 window
            int pixelValues[4];
            pixelValues[0] = input[y * width + x];
            pixelValues[1] = input[y * width + (x + 1)];
            pixelValues[2] = input[(y + 1) * width + x];
            pixelValues[3] = input[(y + 1) * width + (x + 1)];
            
            // Count positive and negative pixel values in the window
            int countPositive = count_if(pixelValues, pixelValues + 4, [](int val) { return val > 0; });
            int countNegative = count_if(pixelValues, pixelValues + 4, [](int val) { return val < 0; });

            // Calculate Laplacian value for the center pixel
            int laplacian = 4 * pixelValues[0] - (pixelValues[1] + pixelValues[2] + pixelValues[3]);

            // If there is a zero-crossing and the Laplacian value is above the threshold, set the pixel to 255, otherwise set it to 0
            output[y * width + x] = (laplacian > threshold) ? 255 : 0;
        }
    }

    // Display a message indicating that Laplacian Operator has been applied and zero-crossings have been found
    cout << "--Laplacian Operator applied and zero-crossings found--" << endl;
} // End of function

// 2. Apply Sobel Operator to the "building.raw" image to obtain the gradient at all pixels
// Function to apply the Sobel operator
void applySobel(const vector<unsigned char>& input, vector<unsigned char>& output, int width, int height) {
    // Sobel kernels for X and Y directions
    int sobelKernelX[3][3] = {{-1, 0, 1},
                               {-2, 0, 2},
                               {-1, 0, 1}};

    int sobelKernelY[3][3] = {{-1, -2, -1},
                               {0, 0, 0},
                               {1, 2, 1}};

    // Ensure the output vector has the correct size
    output.resize(input.size());

    // Apply Sobel operator separately for X and Y directions
    for (int y = 1; y < height - 1; ++y) {
        for (int x = 1; x < width - 1; ++x) {
            // Initialize sums for X and Y gradients
            int sumX = 0;
            int sumY = 0;

            // Iterate over the 3x3 neighborhood of the current pixel
            for (int ky = -1; ky <= 1; ++ky) {
                for (int kx = -1; kx <= 1; ++kx) {
                    // Convolution with Sobel kernels for X and Y directions
                    sumX += sobelKernelX[ky + 1][kx + 1] * input[(y + ky) * width + (x + kx)];
                    sumY += sobelKernelY[ky + 1][kx + 1] * input[(y + ky) * width + (x + kx)];
                }
            }

            // Compute the magnitude of the gradient
            int magnitude = static_cast<int>(sqrt(static_cast<double>(sumX * sumX + sumY * sumY)));

            // Clip the magnitude to 255 (maximum intensity)
            output[y * width + x] = static_cast<unsigned char>(min(255, magnitude));
        }
    }

    // Display a message indicating that Sobel has been applied
    cout << "--Sobel applied--" << endl;
} // End of function

// 3. Generate an edge image where a pixel is an edge point if it is a zero-crossing point 
// And the gradient at this point is greater than or equal to a pre-specified threshold.
// Function to generate the edge image
void generateEdgeImage(const vector<unsigned char>& zeroCrossings, const vector<unsigned char>& gradient, vector<unsigned char>& edgeImage, int width, int height, int threshold) {
    // Ensure the output vector has the correct size
    edgeImage.resize(gradient.size(), 0);

    // Iterate over zero-crossings and apply threshold
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            // Check if gradient at the zero-crossing point exceeds the threshold
            if (zeroCrossings[y * width + x] > 0 && gradient[y * width + x] >= threshold) {
                // Set the pixel in the edge image to a high intensity value
                edgeImage[y * width + x] = 255;
            }
        }
    }

    // Display a message indicating that Edge has been generated
    cout << "--Edge generated--" << endl;
} // End of function

// 4. Implement the Hough transform algorithm to extract three longest linear structures from your edge image.
// Function to implement the Hough transform
void houghTransform(const unsigned char* edgeImage, vector<pair<double, double>>& lines, int width, int height, int threshold) {
    // Calculate the diagonal length of the image
    const int diagonal = static_cast<int>(sqrt(width * width + height * height));
    // Set the maximum theta value for the Hough transform
    const int maxTheta = 180;
    // Set the size of the accumulator array
    const int accumulatorSize = diagonal * 2;
    // Calculate the center coordinates of the image
    const int centerX = width / 2;
    const int centerY = height / 2;

    // Initialize the accumulator array with zeros
    vector<int> accumulator(accumulatorSize * maxTheta, 0);

    // Voting in the accumulator
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            // Check if the pixel is part of an edge
            if (edgeImage[y * width + x] > 0) {
                // Vote for possible lines at different theta values
                for (int theta = 0; theta < maxTheta; ++theta) {
                    int rho = static_cast<int>(x * cos(theta * M_PI / 180.0) + y * sin(theta * M_PI / 180.0)) + diagonal;
                    accumulator[rho * maxTheta + theta]++;
                }
            }
        }
    }

    // Find the peaks in the accumulator and store the corresponding lines
    for (int r = 0; r < accumulatorSize; ++r) {
        for (int theta = 0; theta < maxTheta; ++theta) {
            // Check if the accumulator value is above the threshold
            if (accumulator[r * maxTheta + theta] > threshold) {
                // Store the line parameters (rho and theta) in the lines vector
                lines.emplace_back(r - diagonal, theta);
            }
        }
    }

    // Display a message indicating that Hough Transform is complete
    cout << "--Hough Transform complete--" << endl;
} // End of function

// Function to draw lines on the image
void drawLines(unsigned char* outputImage, vector<pair<double, double>>& lines, int width, int height) {
    // Iterate through each line obtained from the Hough Transform
    for (const auto& line : lines) {
        // Extract rho and theta values
        int rho = line.first;
        int theta = line.second;

        // Draw the line on the output image
        for (int x = 0; x < width; ++x) {
            // Calculate corresponding y coordinate using the line equation
            int y = static_cast<int>((rho - x * cos(theta * M_PI / 180.0)) / sin(theta * M_PI / 180.0));

            // Check if the calculated y coordinate is within the image boundaries
            if (y >= 0 && y < height) {
                // Set the pixel to white (255) to represent the line
                outputImage[y * width + x] = 255;
            }
        }
    }
} // End of function

// Function to save the image data to a RAW file
void saveRawImage(const vector<unsigned char>& imageData, const string& filename, int width, int height) {
    // Open the output file in binary mode
    ofstream file(filename, ios::binary);

    // Check if the file is successfully opened
    if (!file.is_open()) {
        cerr << "Error: Could not open the output file." << endl;
        return;
    }

    // Write the image data to the file
    file.write(reinterpret_cast<const char*>(imageData.data()), imageData.size());

    // Close the file after writing the data
    file.close();
} // End of function


// Main function
int main() {
    // Load the RAW image
    ifstream file("building.raw", ios::binary);

    // Check if the file is successfully opened
    if (!file.is_open()) {
        cerr << "Error: Could not open the input file." << endl;
        return 1;
    }

    // Read the raw image data into a vector
    vector<unsigned char> rawImage((istreambuf_iterator<char>(file)), (istreambuf_iterator<char>()));
    file.close();

    // Replace with the actual width and height of the image
    int width = 560;
    int height = 420;
    int threshold = 110; // Threshold value
    int numLines = 3;    // Number of lines to extract

    // Check if the image dimensions are valid
    if (width <= 0 || height <= 0 || rawImage.size() != width * height) {
        cerr << "Error: Invalid image dimensions or size. Width: " << width << ", Height: " << height << endl;
        return 1;
    }

    // Variables to store intermediate and final results
    vector<unsigned char> laplacianImage;
    vector<pair<int, int>> zeroCrossings;
    vector<unsigned char> gradientImage;
    vector<unsigned char> edgeImage;
    vector<pair<double, double>> lines;

    // Apply Laplacian Operator
    applyLaplacian(rawImage, laplacianImage, width, height, threshold);

    // Save the Laplacian image with zero-crossings
    saveRawImage(laplacianImage, "zero_crossing.raw", width, height);

    // Apply Sobel Operator
    applySobel(rawImage, gradientImage, width, height);

    // Save the gradient image
    saveRawImage(gradientImage, "gradient_image.raw", width, height);

    // Generate Edge Image
    generateEdgeImage(laplacianImage, gradientImage, edgeImage, width, height, threshold);

    // Save the edge image
    saveRawImage(edgeImage, "edge_image.raw", width, height);

    // Hough Transform
    houghTransform(edgeImage.data(), lines, width, height, threshold);

    // Sort lines based on votes (accumulator values)
    sort(lines.begin(), lines.end(), [](const auto& a, const auto& b) {
        return a.second > b.second;
    });

    // Extract the top 3 lines
    lines.resize(min(3, static_cast<int>(lines.size())));

    // Create output image
    vector<unsigned char> outputImage(width * height, 0);
    drawLines(outputImage.data(), lines, width, height);

    // Save the Hough-transformed image with extracted lines
    ofstream outputFile("hough_image.raw", ios::binary);
    if (outputFile.is_open()) {
        outputFile.write(reinterpret_cast<const char*>(outputImage.data()), width * height);
        outputFile.close();
    } else {
        cerr << "Error opening output file." << endl;
        return 1;
    }

    return 0;
} // End of program