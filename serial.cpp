#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <algorithm>
#include <cmath>

const int num_ants = 1;
const double alpha = 1.0;  // Pheromone importance
const double beta = 2.0;   // Heuristic information importance
const double evaporation_rate = 0.5;  // Rate at which pheromone evaporates

class AntColony {
public:
    AntColony(int num_nodes, std::vector<std::vector<double>>& distance_matrix)
        : num_nodes(num_nodes), distance_matrix(distance_matrix) {
        pheromone_matrix = std::vector<std::vector<double>>(num_nodes, std::vector<double>(num_nodes, 1.0));

        // Initialize pheromone_matrix with small random values
        for (auto& row : pheromone_matrix) {
            for (double& val : row) {
                val = (rand() % 10 + 1) / 10.0;
            }
        }
    }

    std::vector<int> find_best_tour() {
        std::vector<int> best_tour;
        double best_distance = std::numeric_limits<double>::infinity();

        for (int i = 0; i < num_ants; ++i) {
            std::vector<int> tour = generate_tour();
            double tour_distance = calculate_tour_distance(tour);

            if (tour_distance < best_distance) {
                best_distance = tour_distance;
                best_tour = tour;
            }

            // Update pheromones based on the tour
            update_pheromones(tour, tour_distance);
        }

        return best_tour;
    }

private:
    int num_nodes;
    std::vector<std::vector<double>> distance_matrix;
    std::vector<std::vector<double>> pheromone_matrix;

    std::vector<int> generate_tour() {
        std::vector<int> tour;
        tour.reserve(num_nodes);

        // Start from a random node
        int current_node = rand() % num_nodes;
        tour.push_back(current_node);

        // Build the tour
        for (int i = 0; i < num_nodes - 1; ++i) {
            current_node = select_next_node(current_node, tour);
            tour.push_back(current_node);
        }

        return tour;
    }

    int select_next_node(int current_node, const std::vector<int>& tour) {
        std::vector<double> probabilities(num_nodes, 0.0);
        double total_prob = 0.0;

        for (int i = 0; i < num_nodes; ++i) {
            if (std::find(tour.begin(), tour.end(), i) == tour.end()) {
                double pheromone = pheromone_matrix[current_node][i];
                double heuristic = 1.0 / distance_matrix[current_node][i];

                probabilities[i] = pow(pheromone, alpha) * pow(heuristic, beta);
                total_prob += probabilities[i];
            }
        }

        // Roulette wheel selection
        double rand_value = (rand() / (RAND_MAX + 1.0)) * total_prob;
        double cumulative_prob = 0.0;

        for (int i = 0; i < num_nodes; ++i) {
            if (std::find(tour.begin(), tour.end(), i) == tour.end()) {
                cumulative_prob += probabilities[i];
                if (cumulative_prob >= rand_value) {
                    return i;
                }
            }
        }

        // Fallback: should not reach here
        return -1;
    }

    double calculate_tour_distance(const std::vector<int>& tour) {
        double distance = 0.0;
        for (int i = 0; i < num_nodes - 1; ++i) {
            distance += distance_matrix[tour[i]][tour[i + 1]];
        }
        // Return to the starting node
        distance += distance_matrix[tour.back()][tour.front()];
        return distance;
    }

    void update_pheromones(const std::vector<int>& tour, double tour_distance) {
        for (int i = 0; i < num_nodes - 1; ++i) {
            int from = tour[i];
            int to = tour[i + 1];
            pheromone_matrix[from][to] = (1.0 - evaporation_rate) * pheromone_matrix[from][to] + (1.0 / tour_distance);
            pheromone_matrix[to][from] = pheromone_matrix[from][to];  // Symmetric TSP
        }
    }
};

int main() {
    // Create a 100x100 distance matrix with random values
    srand(static_cast<unsigned>(time(nullptr)));
    const int num_nodes = 1000;
    std::vector<std::vector<double>> distance_matrix(num_nodes, std::vector<double>(num_nodes, 0.0));
    for (int i = 0; i < num_nodes; ++i) {
        for (int j = 0; j < num_nodes; ++j) {
            if (i != j) {
                distance_matrix[i][j] = rand() % 100 + 1;  // Random distance between 1 and 100
            }
        }
    }

    AntColony ant_colony(num_nodes, distance_matrix);

    double start_time = clock();
    std::vector<int> best_tour = ant_colony.find_best_tour();
    double end_time = clock();

    // Print the best tour
    std::cout << "Best Tour: ";
    for (int node : best_tour) {
        std::cout << node << " ";
    }
    std::cout << "\n";

    // Print the runtime
    std::cout << "Runtime: " << end_time - start_time << " milliseconds\n";

    return 0;
}
