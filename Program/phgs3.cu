//
// Created by chkwon on 3/23/22.
//

#include "AlgorithmParameters.h"
#include <iostream>

extern "C" struct AlgorithmParameters default_algorithm_parameters()
{
    struct AlgorithmParameters ap
    {
    };

    ap.nbGranular = 20;
    ap.mu = 25;
    ap.lambda = 40;
    ap.nbElite = 4;
    ap.nbClose = 5;

    ap.nbIterPenaltyManagement = 100;
    ap.targetFeasible = 0.2;
    ap.penaltyDecrease = 0.85;
    ap.penaltyIncrease = 1.2;

    ap.seed = 0;
    ap.nbIter = 20000;
    ap.nbIterTraces = 500;
    ap.timeLimit = 0;
    ap.useSwapStar = 1;

    return ap;
}

void print_algorithm_parameters(const AlgorithmParameters &ap)
{
    std::cout << "=========== Algorithm Parameters =================" << std::endl;
    std::cout << "---- nbGranular              is set to " << ap.nbGranular << std::endl;
    std::cout << "---- mu                      is set to " << ap.mu << std::endl;
    std::cout << "---- lambda                  is set to " << ap.lambda << std::endl;
    std::cout << "---- nbElite                 is set to " << ap.nbElite << std::endl;
    std::cout << "---- nbClose                 is set to " << ap.nbClose << std::endl;
    std::cout << "---- nbIterPenaltyManagement is set to " << ap.nbIterPenaltyManagement << std::endl;
    std::cout << "---- targetFeasible          is set to " << ap.targetFeasible << std::endl;
    std::cout << "---- penaltyDecrease         is set to " << ap.penaltyDecrease << std::endl;
    std::cout << "---- penaltyIncrease         is set to " << ap.penaltyIncrease << std::endl;
    std::cout << "---- seed                    is set to " << ap.seed << std::endl;
    std::cout << "---- nbIter                  is set to " << ap.nbIter << std::endl;
    std::cout << "---- nbIterTraces            is set to " << ap.nbIterTraces << std::endl;
    std::cout << "---- timeLimit               is set to " << ap.timeLimit << std::endl;
    std::cout << "---- useSwapStar             is set to " << ap.useSwapStar << std::endl;
    std::cout << "==================================================" << std::endl;
}

//
// Created by chkwon on 3/23/22.
//

#include "C_Interface.h"
#include "Population.h"
#include "Params.h"
#include "Genetic.h"
#include <string>
#include <iostream>
#include <vector>
#include <cmath>

Solution *prepare_solution(Population &population, Params &params)
{
    // Preparing the best solution
    Solution *sol = new Solution;
    sol->time = (double)(clock() - params.startTime) / (double)CLOCKS_PER_SEC;

    if (population.getBestFound() != nullptr)
    {
        // Best individual
        auto best = population.getBestFound();

        // setting the cost
        sol->cost = best->eval.penalizedCost;

        // finding out the number of routes in the best individual
        int n_routes = 0;
        // bookmarkreduction
        for (int k = 0; k < params.nbVehicles; k++)
            if (!best->chromR[k].empty())
                ++n_routes;

        // filling out the route information
        sol->n_routes = n_routes;
        sol->routes = new SolutionRoute[n_routes];
        // bookmarkimp
        for (int k = 0; k < n_routes; k++)
        {
            sol->routes[k].length = (int)best->chromR[k].size();
            sol->routes[k].path = new int[sol->routes[k].length];
            std::copy(best->chromR[k].begin(), best->chromR[k].end(), sol->routes[k].path);
        }
    }
    else
    {
        sol->cost = 0.0;
        sol->n_routes = 0;
        sol->routes = nullptr;
    }
    return sol;
}

extern "C" Solution *solve_cvrp(
    int n, double *x, double *y, double *serv_time, double *dem,
    double vehicleCapacity, double durationLimit, char isRoundingInteger, char isDurationConstraint,
    int max_nbVeh, const AlgorithmParameters *ap, char verbose)
{
    Solution *result;

    try
    {
        std::vector<double> x_coords(x, x + n);
        std::vector<double> y_coords(y, y + n);
        std::vector<double> service_time(serv_time, serv_time + n);
        std::vector<double> demands(dem, dem + n);

        std::vector<std::vector<double>> distance_matrix(n, std::vector<double>(n));
        // bookmarkimp
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                distance_matrix[i][j] = std::sqrt(
                    (x_coords[i] - x_coords[j]) * (x_coords[i] - x_coords[j]) + (y_coords[i] - y_coords[j]) * (y_coords[i] - y_coords[j]));
                if (isRoundingInteger)
                    distance_matrix[i][j] = std::round(distance_matrix[i][j]);
            }
        }

        Params params(x_coords, y_coords, distance_matrix, service_time, demands, vehicleCapacity, durationLimit, max_nbVeh, isDurationConstraint, verbose, *ap);

        // Running HGS and returning the result
        Genetic solver(params);
        solver.run();
        result = prepare_solution(solver.population, params);
    }
    catch (const std::string &e)
    {
        std::cout << "EXCEPTION | " << e << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cout << "EXCEPTION | " << e.what() << std::endl;
    }

    return result;
}

extern "C" Solution *solve_cvrp_dist_mtx(
    int n, double *x, double *y, double *dist_mtx, double *serv_time, double *dem,
    double vehicleCapacity, double durationLimit, char isDurationConstraint,
    int max_nbVeh, const AlgorithmParameters *ap, char verbose)
{
    Solution *result;
    std::vector<double> x_coords;
    std::vector<double> y_coords;

    try
    {
        if (x != nullptr && y != nullptr)
        {
            x_coords = {x, x + n};
            y_coords = {y, y + n};
        }

        std::vector<double> service_time(serv_time, serv_time + n);
        std::vector<double> demands(dem, dem + n);

        std::vector<std::vector<double>> distance_matrix(n, std::vector<double>(n));
        for (int i = 0; i < n; i++)
        { // row
            for (int j = 0; j < n; j++)
            { // column
                distance_matrix[i][j] = dist_mtx[n * i + j];
            }
        }

        Params params(x_coords, y_coords, distance_matrix, service_time, demands, vehicleCapacity, durationLimit, max_nbVeh, isDurationConstraint, verbose, *ap);

        // Running HGS and returning the result
        Genetic solver(params);
        solver.run();
        result = prepare_solution(solver.population, params);
    }
    catch (const std::string &e)
    {
        std::cout << "EXCEPTION | " << e << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cout << "EXCEPTION | " << e.what() << std::endl;
    }

    return result;
}

extern "C" void delete_solution(Solution *sol)
{
    for (int i = 0; i < sol->n_routes; ++i)
        delete[] sol->routes[i].path;

    delete[] sol->routes;
    delete sol;
}

#include "Genetic.h"
#define NUM_THREADS 1024
#define BLOCKS 1024

void Genetic::run()
{
    /* INITIAL POPULATION */
    population.generatePopulation();

    int nbIter;
    int nbIterNonProd = 1;
    if (params.verbose)
        std::cout << "----- STARTING GENETIC ALGORITHM" << std::endl;
    for (nbIter = 0; nbIterNonProd <= params.ap.nbIter && (params.ap.timeLimit == 0 || (double)(clock() - params.startTime) < params.ap.timeLimit); nbIter++)
    {
        /* SELECTION AND CROSSOVER */
        crossoverOX(offspring, population.getBinaryTournament(), population.getBinaryTournament());

        /*[edit] LOCAL SEARCH */
        localSearch.run(offspring, params.penaltyCapacity, params.penaltyDuration);

        bool isNewBest = population.addIndividual(offspring, true);
        if (!offspring.eval.isFeasible && params.ran() % 2 == 0) // Repair half of the solutions in case of infeasibility
        {
            localSearch.run(offspring, params.penaltyCapacity * 10., params.penaltyDuration * 10.);
            if (offspring.eval.isFeasible)
                isNewBest = (population.addIndividual(offspring, false) || isNewBest);
        }

        /* TRACKING THE NUMBER OF ITERATIONS SINCE LAST SOLUTION IMPROVEMENT */
        if (isNewBest)
            nbIterNonProd = 1;
        else
            nbIterNonProd++;

        /* DIVERSIFICATION, PENALTY MANAGEMENT AND TRACES */
        if (nbIter % params.ap.nbIterPenaltyManagement == 0)
            population.managePenalties();
        if (nbIter % params.ap.nbIterTraces == 0)
            population.printState(nbIter, nbIterNonProd);

        /* FOR TESTS INVOLVING SUCCESSIVE RUNS UNTIL A TIME LIMIT: WE RESET THE ALGORITHM/POPULATION EACH TIME maxIterNonProd IS ATTAINED*/
        if (params.ap.timeLimit != 0 && nbIterNonProd == params.ap.nbIter)
        {
            population.restart();
            nbIterNonProd = 1;
        }
    }
    if (params.verbose)
        std::cout << "----- GENETIC ALGORITHM FINISHED AFTER " << nbIter << " ITERATIONS. TIME SPENT: " << (double)(clock() - params.startTime) / (double)CLOCKS_PER_SEC << std::endl;
}

__global__ void crossoverOX_kernel(int nbClients, bool *freqClient, int *parent2ChromT, int end, int *resultChromT, int j)
{
    for (int i = 1; i <= nbClients; i++)
    {
        int temp = parent2ChromT[(end + i) % nbClients];
        if (freqClient[temp] == false)
        {
            resultChromT[j % nbClients] = temp;
            j++;
        }
    }
}

void Genetic::crossoverOX(Individual &result, const Individual &parent1, const Individual &parent2)
{
    // Frequency table to track the customers which have been already inserted
    std::vector<bool> freqClient = std::vector<bool>(params.nbClients + 1, false);

    // Picking the beginning and end of the crossover zone
    std::uniform_int_distribution<> distr(0, params.nbClients - 1);
    int start = distr(params.ran);
    int end = distr(params.ran);

    // Avoid that start and end coincide by accident
    while (end == start)
        end = distr(params.ran);

    // Copy from start to end
    int j = start, j2 = start;
    while (j % params.nbClients != (end + 1) % params.nbClients)
    {
        result.chromT[j % params.nbClients] = parent1.chromT[j % params.nbClients];
        freqClient[result.chromT[j % params.nbClients]] = true;
        j++;
        j2++;
    }
    //[bookmark]
    // Fill the remaining elements in the order given by the second parent

    bool *freqClient2 = (bool *)malloc(sizeof(bool) * params.nbClients);
    int *parent2ChromT2 = (int *)malloc(sizeof(int) * params.nbClients);
    int *resultChromT2 = (int *)malloc(sizeof(int) * params.nbClients);

    bool *parallel_freqClient;
    int *parallel_parent2ChromT;
    int *parallel_resultChromT;

    int count = 0;
    for (int i = 1; i <= params.nbClients; i++)
    {
        count++;
        parent2ChromT2[(end + i) % params.nbClients] = parent2.chromT[(end + i) % params.nbClients];
        freqClient2[parent2.chromT[(end + i) % params.nbClients]] = freqClient[parent2.chromT[(end + i) % params.nbClients]];
        {
            resultChromT2[j % params.nbClients] = result.chromT[j % params.nbClients];
            j++;
        }
    }

    j = j2;

    cudaMalloc((void **)&parallel_freqClient, count * sizeof(bool));
    cudaMalloc((void **)&parallel_parent2ChromT, count * sizeof(Individual));
    cudaMalloc((void **)&parallel_resultChromT, count * sizeof(Individual));

    cudaMemcpy(parallel_freqClient, freqClient2, count * sizeof(bool), cudaMemcpyHostToDevice);
    cudaMemcpy(parallel_parent2ChromT, parent2ChromT2, count * sizeof(Individual), cudaMemcpyHostToDevice);
    cudaMemcpy(parallel_resultChromT, resultChromT2, count * sizeof(Individual), cudaMemcpyHostToDevice);

    int divisions = (params.nbClients / BLOCKS) / NUM_THREADS;
    crossoverOX_kernel<<<BLOCKS, NUM_THREADS>>>(params.nbClients, parallel_freqClient, parallel_parent2ChromT, end, parallel_resultChromT, j);

    cudaMemcpy(parallel_freqClient, freqClient2, count * sizeof(bool), cudaMemcpyDeviceToHost);
    cudaMemcpy(parallel_parent2ChromT, parent2ChromT2, count * sizeof(Individual), cudaMemcpyDeviceToHost);
    cudaMemcpy(parallel_resultChromT, resultChromT2, count * sizeof(Individual), cudaMemcpyDeviceToHost);

    j = j2;

    for (int i = 1; i <= params.nbClients; i++)
    {
        int temp = parent2.chromT[(end + i) % params.nbClients];
        if (freqClient[temp] == false)
        {
            result.chromT[j % params.nbClients] = resultChromT2[j % params.nbClients];
            j++;
        }
    }

    // for (int i = 1; i <= params.nbClients; i++)
    // {
    // int temp = parent2.chromT[(end + i) % params.nbClients];
    // if (freqClient[temp] == false)
    // {
    //     result.chromT[j % params.nbClients] = temp;
    //     j++;
    // }
    // }

    // Complete the individual with the Split algorithm
    split.generalSplit(result, parent1.eval.nbRoutes);
}

Genetic::Genetic(Params &params) : params(params),
                                   split(params),
                                   localSearch(params),
                                   population(params, this->split, this->localSearch),
                                   offspring(params) {}

#include "Individual.h"

void Individual::evaluateCompleteCost(const Params &params)
{
    eval = EvalIndiv();
    // bookmarkimpbut has paramas.timeCost
    for (int r = 0; r < params.nbVehicles; r++)
    {
        if (!chromR[r].empty())
        {
            double distance = params.timeCost[0][chromR[r][0]];
            double load = params.cli[chromR[r][0]].demand;
            double service = params.cli[chromR[r][0]].serviceDuration;
            predecessors[chromR[r][0]] = 0;
            for (int i = 1; i < (int)chromR[r].size(); i++)
            {
                distance += params.timeCost[chromR[r][i - 1]][chromR[r][i]];
                load += params.cli[chromR[r][i]].demand;
                service += params.cli[chromR[r][i]].serviceDuration;
                predecessors[chromR[r][i]] = chromR[r][i - 1];
                successors[chromR[r][i - 1]] = chromR[r][i];
            }
            successors[chromR[r][chromR[r].size() - 1]] = 0;
            distance += params.timeCost[chromR[r][chromR[r].size() - 1]][0];
            eval.distance += distance;
            eval.nbRoutes++;
            if (load > params.vehicleCapacity)
                eval.capacityExcess += load - params.vehicleCapacity;
            if (distance + service > params.durationLimit)
                eval.durationExcess += distance + service - params.durationLimit;
        }
    }

    eval.penalizedCost = eval.distance + eval.capacityExcess * params.penaltyCapacity + eval.durationExcess * params.penaltyDuration;
    eval.isFeasible = (eval.capacityExcess < MY_EPSILON && eval.durationExcess < MY_EPSILON);
}

Individual::Individual(Params &params)
{
    successors = std::vector<int>(params.nbClients + 1);
    predecessors = std::vector<int>(params.nbClients + 1);
    chromR = std::vector<std::vector<int>>(params.nbVehicles);
    chromT = std::vector<int>(params.nbClients);
    // bookmark
    for (int i = 0; i < params.nbClients; i++)
        chromT[i] = i + 1;
    std::shuffle(chromT.begin(), chromT.end(), params.ran);
    eval.penalizedCost = 1.e30;
}

Individual::Individual(Params &params, std::string fileName) : Individual(params)
{
    double readCost;
    chromT.clear();
    std::ifstream inputFile(fileName);
    if (inputFile.is_open())
    {
        std::string inputString;
        inputFile >> inputString;
        // Loops in the input file as long as the first line keyword is "Route"
        for (int r = 0; inputString == "Route"; r++)
        {
            inputFile >> inputString;
            getline(inputFile, inputString);
            std::stringstream ss(inputString);
            int inputCustomer;
            while (ss >> inputCustomer) // Loops as long as there is an integer to read in this route
            {
                chromT.push_back(inputCustomer);
                chromR[r].push_back(inputCustomer);
            }
            inputFile >> inputString;
        }
        if (inputString == "Cost")
            inputFile >> readCost;
        else
            throw std::string("Unexpected token in input solution");

        // Some safety checks and printouts
        evaluateCompleteCost(params);
        if ((int)chromT.size() != params.nbClients)
            throw std::string("Input solution does not contain the correct number of clients");
        if (!eval.isFeasible)
            throw std::string("Input solution is infeasible");
        if (eval.penalizedCost != readCost)
            throw std::string("Input solution has a different cost than announced in the file");
        if (params.verbose)
            std::cout << "----- INPUT SOLUTION HAS BEEN SUCCESSFULLY READ WITH COST " << eval.penalizedCost << std::endl;
    }
    else
        throw std::string("Impossible to open solution file provided in input in : " + fileName);
}

//
// Created by chkwon on 3/22/22.
//

#include <fstream>
#include <cmath>
#include "InstanceCVRPLIB.h"

InstanceCVRPLIB::InstanceCVRPLIB(std::string pathToInstance, bool isRoundingInteger = true)
{
    std::string content, content2, content3;
    double serviceTimeData = 0.;

    // Read INPUT dataset
    std::ifstream inputFile(pathToInstance);
    if (inputFile.is_open())
    {
        getline(inputFile, content);
        getline(inputFile, content);
        getline(inputFile, content);
        for (inputFile >> content; content != "NODE_COORD_SECTION"; inputFile >> content)
        {
            if (content == "DIMENSION")
            {
                inputFile >> content2 >> nbClients;
                nbClients--;
            } // Need to substract the depot from the number of nodes
            else if (content == "EDGE_WEIGHT_TYPE")
                inputFile >> content2 >> content3;
            else if (content == "CAPACITY")
                inputFile >> content2 >> vehicleCapacity;
            else if (content == "DISTANCE")
            {
                inputFile >> content2 >> durationLimit;
                isDurationConstraint = true;
            }
            else if (content == "SERVICE_TIME")
                inputFile >> content2 >> serviceTimeData;
            else
                throw std::string("Unexpected data in input file: " + content);
        }
        if (nbClients <= 0)
            throw std::string("Number of nodes is undefined");
        if (vehicleCapacity == 1.e30)
            throw std::string("Vehicle capacity is undefined");

        x_coords = std::vector<double>(nbClients + 1);
        y_coords = std::vector<double>(nbClients + 1);
        demands = std::vector<double>(nbClients + 1);
        service_time = std::vector<double>(nbClients + 1);

        // Reading node coordinates
        // depot must be the first element
        // 		- i = 0 in the for-loop below, or
        // 		- node_number = 1 in the .vrp file
        // customers are
        // 		- i = 1, 2, ..., nbClients in the for-loop below, or
        // 		- node_number = 2, 3, ..., nb_Clients in the .vrp file
        int node_number;
        for (int i = 0; i <= nbClients; i++)
        {
            inputFile >> node_number >> x_coords[i] >> y_coords[i];
            if (node_number != i + 1)
                throw std::string("The node numbering is not in order.");
        }

        // Reading demand information
        inputFile >> content;
        if (content != "DEMAND_SECTION")
            throw std::string("Unexpected data in input file: " + content);
        for (int i = 0; i <= nbClients; i++)
        {
            inputFile >> content >> demands[i];
            service_time[i] = (i == 0) ? 0. : serviceTimeData;
        }

        // Calculating 2D Euclidean Distance
        dist_mtx = std::vector<std::vector<double>>(nbClients + 1, std::vector<double>(nbClients + 1));
        for (int i = 0; i <= nbClients; i++)
        {
            for (int j = 0; j <= nbClients; j++)
            {
                dist_mtx[i][j] = std::sqrt(
                    (x_coords[i] - x_coords[j]) * (x_coords[i] - x_coords[j]) + (y_coords[i] - y_coords[j]) * (y_coords[i] - y_coords[j]));

                if (isRoundingInteger)
                    dist_mtx[i][j] = round(dist_mtx[i][j]);
            }
        }

        // Reading depot information (in all current instances the depot is represented as node 1, the program will return an error otherwise)
        inputFile >> content >> content2 >> content3 >> content3;
        if (content != "DEPOT_SECTION")
            throw std::string("Unexpected data in input file: " + content);
        if (content2 != "1")
            throw std::string("Expected depot index 1 instead of " + content2);
        if (content3 != "EOF")
            throw std::string("Unexpected data in input file: " + content3);
    }
    else
        throw std::string("Impossible to open instance file: " + pathToInstance);
}

#include "LocalSearch.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#define NUM_THREADS 1024
#define BLOCKS 1024

void LocalSearch::run(Individual &indiv, double penaltyCapacityLS, double penaltyDurationLS)
{
    this->penaltyCapacityLS = penaltyCapacityLS;
    this->penaltyDurationLS = penaltyDurationLS;
    loadIndividual(indiv);

    // Shuffling the order of the nodes explored by the LS to allow for more diversity in the search
    std::shuffle(orderNodes.begin(), orderNodes.end(), params.ran);
    std::shuffle(orderRoutes.begin(), orderRoutes.end(), params.ran);

    for (int i = 1; i <= params.nbClients; i++)
        if (params.ran() % params.ap.nbGranular == 0) // O(n/nbGranular) calls to the inner function on average, to achieve linear-time complexity overall
            std::shuffle(params.correlatedVertices[i].begin(), params.correlatedVertices[i].end(), params.ran);

    searchCompleted = false;
    for (loopID = 0; !searchCompleted; loopID++)
    {
        if (loopID > 1) // Allows at least two loops since some moves involving empty routes are not checked at the first loop
            searchCompleted = true;

        //[edit]parallelizing from here

        /*[edit] CLASSICAL ROUTE IMPROVEMENT (RI) MOVES SUBJECT TO A PROXIMITY RESTRICTION */
        for (int posU = 0; posU < params.nbClients; posU++)
        {
            nodeU = &clients[orderNodes[posU]]; //[edit]
            int lastTestRINodeU = nodeU->whenLastTestedRI;
            nodeU->whenLastTestedRI = nbMoves;
            for (int posV = 0; posV < (int)params.correlatedVertices[nodeU->cour].size(); posV++)
            {
                nodeV = &clients[params.correlatedVertices[nodeU->cour][posV]];
                if (loopID == 0 || std::max<int>(nodeU->route->whenLastModified, nodeV->route->whenLastModified) > lastTestRINodeU) // only evaluate moves involving routes that have been modified since last move evaluations for nodeU
                {
                    // Randomizing the order of the neighborhoods within this loop does not matter much as we are already randomizing the order of the node pairs (and it's not very common to find improving moves of different types for the same node pair)
                    setLocalVariablesRouteU();
                    setLocalVariablesRouteV();
                    if (move1())
                        continue; // RELOCATE
                    if (move2())
                        continue; // RELOCATE
                    if (move3())
                        continue; // RELOCATE
                    if (nodeUIndex <= nodeVIndex && move4())
                        continue; // SWAP
                    if (move5())
                        continue; // SWAP
                    if (nodeUIndex <= nodeVIndex && move6())
                        continue; // SWAP
                    if (intraRouteMove && move7())
                        continue; // 2-OPT
                    if (!intraRouteMove && move8())
                        continue; // 2-OPT*
                    if (!intraRouteMove && move9())
                        continue; // 2-OPT*

                    // Trying moves that insert nodeU directly after the depot
                    if (nodeV->prev->isDepot)
                    {
                        nodeV = nodeV->prev;
                        setLocalVariablesRouteV();
                        if (move1())
                            continue; // RELOCATE
                        if (move2())
                            continue; // RELOCATE
                        if (move3())
                            continue; // RELOCATE
                        if (!intraRouteMove && move8())
                            continue; // 2-OPT*
                        if (!intraRouteMove && move9())
                            continue; // 2-OPT*
                    }
                }
            }

            /* MOVES INVOLVING AN EMPTY ROUTE -- NOT TESTED IN THE FIRST LOOP TO AVOID INCREASING TOO MUCH THE FLEET SIZE */
            if (loopID > 0 && !emptyRoutes.empty())
            {
                nodeV = routes[*emptyRoutes.begin()].depot;
                setLocalVariablesRouteU();
                setLocalVariablesRouteV();
                if (move1())
                    continue; // RELOCATE
                if (move2())
                    continue; // RELOCATE
                if (move3())
                    continue; // RELOCATE
                if (move9())
                    continue; // 2-OPT*
            }
        }
    }
    if (params.ap.useSwapStar == 1 && params.areCoordinatesProvided)
    {
        /* (SWAP*) MOVES LIMITED TO ROUTE PAIRS WHOSE CIRCLE SECTORS OVERLAP */
        for (int rU = 0; rU < params.nbVehicles; rU++)
        {
            routeU = &routes[orderRoutes[rU]];
            int lastTestSWAPStarRouteU = routeU->whenLastTestedSWAPStar;
            routeU->whenLastTestedSWAPStar = nbMoves;
            for (int rV = 0; rV < params.nbVehicles; rV++)
            {
                routeV = &routes[orderRoutes[rV]];
                if (routeU->nbCustomers > 0 && routeV->nbCustomers > 0 && routeU->cour < routeV->cour && (loopID == 0 || std::max<int>(routeU->whenLastModified, routeV->whenLastModified) > lastTestSWAPStarRouteU))
                    if (CircleSector::overlap(routeU->sector, routeV->sector))
                        swapStar();
            }
        }
    }

    // Register the solution produced by the LS in the individual
    exportIndividual(indiv);
}

void LocalSearch::setLocalVariablesRouteU()
{
    routeU = nodeU->route;
    nodeX = nodeU->next;
    nodeXNextIndex = nodeX->next->cour;
    nodeUIndex = nodeU->cour;
    nodeUPrevIndex = nodeU->prev->cour;
    nodeXIndex = nodeX->cour;
    loadU = params.cli[nodeUIndex].demand;
    serviceU = params.cli[nodeUIndex].serviceDuration;
    loadX = params.cli[nodeXIndex].demand;
    serviceX = params.cli[nodeXIndex].serviceDuration;
}

void LocalSearch::setLocalVariablesRouteV()
{
    routeV = nodeV->route;
    nodeY = nodeV->next;
    nodeYNextIndex = nodeY->next->cour;
    nodeVIndex = nodeV->cour;
    nodeVPrevIndex = nodeV->prev->cour;
    nodeYIndex = nodeY->cour;
    loadV = params.cli[nodeVIndex].demand;
    serviceV = params.cli[nodeVIndex].serviceDuration;
    loadY = params.cli[nodeYIndex].demand;
    serviceY = params.cli[nodeYIndex].serviceDuration;
    intraRouteMove = (routeU == routeV);
}

bool LocalSearch::move1()
{
    double costSuppU = params.timeCost[nodeUPrevIndex][nodeXIndex] - params.timeCost[nodeUPrevIndex][nodeUIndex] - params.timeCost[nodeUIndex][nodeXIndex];
    double costSuppV = params.timeCost[nodeVIndex][nodeUIndex] + params.timeCost[nodeUIndex][nodeYIndex] - params.timeCost[nodeVIndex][nodeYIndex];

    if (!intraRouteMove)
    {
        // Early move pruning to save CPU time. Guarantees that this move cannot improve without checking additional (load, duration...) constraints
        if (costSuppU + costSuppV >= routeU->penalty + routeV->penalty)
            return false;

        costSuppU += penaltyExcessDuration(routeU->duration + costSuppU - serviceU) + penaltyExcessLoad(routeU->load - loadU) - routeU->penalty;

        costSuppV += penaltyExcessDuration(routeV->duration + costSuppV + serviceU) + penaltyExcessLoad(routeV->load + loadU) - routeV->penalty;
    }

    if (costSuppU + costSuppV > -MY_EPSILON)
        return false;
    if (nodeUIndex == nodeYIndex)
        return false;

    insertNode(nodeU, nodeV);
    nbMoves++; // Increment move counter before updating route data
    searchCompleted = false;
    updateRouteData(routeU);
    if (!intraRouteMove)
        updateRouteData(routeV);
    return true;
}

bool LocalSearch::move2()
{
    double costSuppU = params.timeCost[nodeUPrevIndex][nodeXNextIndex] - params.timeCost[nodeUPrevIndex][nodeUIndex] - params.timeCost[nodeXIndex][nodeXNextIndex];
    double costSuppV = params.timeCost[nodeVIndex][nodeUIndex] + params.timeCost[nodeXIndex][nodeYIndex] - params.timeCost[nodeVIndex][nodeYIndex];

    if (!intraRouteMove)
    {
        // Early move pruning to save CPU time. Guarantees that this move cannot improve without checking additional (load, duration...) constraints
        if (costSuppU + costSuppV >= routeU->penalty + routeV->penalty)
            return false;

        costSuppU += penaltyExcessDuration(routeU->duration + costSuppU - params.timeCost[nodeUIndex][nodeXIndex] - serviceU - serviceX) + penaltyExcessLoad(routeU->load - loadU - loadX) - routeU->penalty;

        costSuppV += penaltyExcessDuration(routeV->duration + costSuppV + params.timeCost[nodeUIndex][nodeXIndex] + serviceU + serviceX) + penaltyExcessLoad(routeV->load + loadU + loadX) - routeV->penalty;
    }

    if (costSuppU + costSuppV > -MY_EPSILON)
        return false;
    if (nodeU == nodeY || nodeV == nodeX || nodeX->isDepot)
        return false;

    insertNode(nodeU, nodeV);
    insertNode(nodeX, nodeU);
    nbMoves++; // Increment move counter before updating route data
    searchCompleted = false;
    updateRouteData(routeU);
    if (!intraRouteMove)
        updateRouteData(routeV);
    return true;
}

bool LocalSearch::move3()
{
    double costSuppU = params.timeCost[nodeUPrevIndex][nodeXNextIndex] - params.timeCost[nodeUPrevIndex][nodeUIndex] - params.timeCost[nodeUIndex][nodeXIndex] - params.timeCost[nodeXIndex][nodeXNextIndex];
    double costSuppV = params.timeCost[nodeVIndex][nodeXIndex] + params.timeCost[nodeXIndex][nodeUIndex] + params.timeCost[nodeUIndex][nodeYIndex] - params.timeCost[nodeVIndex][nodeYIndex];

    if (!intraRouteMove)
    {
        // Early move pruning to save CPU time. Guarantees that this move cannot improve without checking additional (load, duration...) constraints
        if (costSuppU + costSuppV >= routeU->penalty + routeV->penalty)
            return false;

        costSuppU += penaltyExcessDuration(routeU->duration + costSuppU - serviceU - serviceX) + penaltyExcessLoad(routeU->load - loadU - loadX) - routeU->penalty;

        costSuppV += penaltyExcessDuration(routeV->duration + costSuppV + serviceU + serviceX) + penaltyExcessLoad(routeV->load + loadU + loadX) - routeV->penalty;
    }

    if (costSuppU + costSuppV > -MY_EPSILON)
        return false;
    if (nodeU == nodeY || nodeX == nodeV || nodeX->isDepot)
        return false;

    insertNode(nodeX, nodeV);
    insertNode(nodeU, nodeX);
    nbMoves++; // Increment move counter before updating route data
    searchCompleted = false;
    updateRouteData(routeU);
    if (!intraRouteMove)
        updateRouteData(routeV);
    return true;
}

bool LocalSearch::move4()
{
    double costSuppU = params.timeCost[nodeUPrevIndex][nodeVIndex] + params.timeCost[nodeVIndex][nodeXIndex] - params.timeCost[nodeUPrevIndex][nodeUIndex] - params.timeCost[nodeUIndex][nodeXIndex];
    double costSuppV = params.timeCost[nodeVPrevIndex][nodeUIndex] + params.timeCost[nodeUIndex][nodeYIndex] - params.timeCost[nodeVPrevIndex][nodeVIndex] - params.timeCost[nodeVIndex][nodeYIndex];

    if (!intraRouteMove)
    {
        // Early move pruning to save CPU time. Guarantees that this move cannot improve without checking additional (load, duration...) constraints
        if (costSuppU + costSuppV >= routeU->penalty + routeV->penalty)
            return false;

        costSuppU += penaltyExcessDuration(routeU->duration + costSuppU + serviceV - serviceU) + penaltyExcessLoad(routeU->load + loadV - loadU) - routeU->penalty;

        costSuppV += penaltyExcessDuration(routeV->duration + costSuppV - serviceV + serviceU) + penaltyExcessLoad(routeV->load + loadU - loadV) - routeV->penalty;
    }

    if (costSuppU + costSuppV > -MY_EPSILON)
        return false;
    if (nodeUIndex == nodeVPrevIndex || nodeUIndex == nodeYIndex)
        return false;

    swapNode(nodeU, nodeV);
    nbMoves++; // Increment move counter before updating route data
    searchCompleted = false;
    updateRouteData(routeU);
    if (!intraRouteMove)
        updateRouteData(routeV);
    return true;
}

bool LocalSearch::move5()
{
    double costSuppU = params.timeCost[nodeUPrevIndex][nodeVIndex] + params.timeCost[nodeVIndex][nodeXNextIndex] - params.timeCost[nodeUPrevIndex][nodeUIndex] - params.timeCost[nodeXIndex][nodeXNextIndex];
    double costSuppV = params.timeCost[nodeVPrevIndex][nodeUIndex] + params.timeCost[nodeXIndex][nodeYIndex] - params.timeCost[nodeVPrevIndex][nodeVIndex] - params.timeCost[nodeVIndex][nodeYIndex];

    if (!intraRouteMove)
    {
        // Early move pruning to save CPU time. Guarantees that this move cannot improve without checking additional (load, duration...) constraints
        if (costSuppU + costSuppV >= routeU->penalty + routeV->penalty)
            return false;

        costSuppU += penaltyExcessDuration(routeU->duration + costSuppU - params.timeCost[nodeUIndex][nodeXIndex] + serviceV - serviceU - serviceX) + penaltyExcessLoad(routeU->load + loadV - loadU - loadX) - routeU->penalty;

        costSuppV += penaltyExcessDuration(routeV->duration + costSuppV + params.timeCost[nodeUIndex][nodeXIndex] - serviceV + serviceU + serviceX) + penaltyExcessLoad(routeV->load + loadU + loadX - loadV) - routeV->penalty;
    }

    if (costSuppU + costSuppV > -MY_EPSILON)
        return false;
    if (nodeU == nodeV->prev || nodeX == nodeV->prev || nodeU == nodeY || nodeX->isDepot)
        return false;

    swapNode(nodeU, nodeV);
    insertNode(nodeX, nodeU);
    nbMoves++; // Increment move counter before updating route data
    searchCompleted = false;
    updateRouteData(routeU);
    if (!intraRouteMove)
        updateRouteData(routeV);
    return true;
}

bool LocalSearch::move6()
{
    double costSuppU = params.timeCost[nodeUPrevIndex][nodeVIndex] + params.timeCost[nodeYIndex][nodeXNextIndex] - params.timeCost[nodeUPrevIndex][nodeUIndex] - params.timeCost[nodeXIndex][nodeXNextIndex];
    double costSuppV = params.timeCost[nodeVPrevIndex][nodeUIndex] + params.timeCost[nodeXIndex][nodeYNextIndex] - params.timeCost[nodeVPrevIndex][nodeVIndex] - params.timeCost[nodeYIndex][nodeYNextIndex];

    if (!intraRouteMove)
    {
        // Early move pruning to save CPU time. Guarantees that this move cannot improve without checking additional (load, duration...) constraints
        if (costSuppU + costSuppV >= routeU->penalty + routeV->penalty)
            return false;

        costSuppU += penaltyExcessDuration(routeU->duration + costSuppU - params.timeCost[nodeUIndex][nodeXIndex] + params.timeCost[nodeVIndex][nodeYIndex] + serviceV + serviceY - serviceU - serviceX) + penaltyExcessLoad(routeU->load + loadV + loadY - loadU - loadX) - routeU->penalty;

        costSuppV += penaltyExcessDuration(routeV->duration + costSuppV + params.timeCost[nodeUIndex][nodeXIndex] - params.timeCost[nodeVIndex][nodeYIndex] - serviceV - serviceY + serviceU + serviceX) + penaltyExcessLoad(routeV->load + loadU + loadX - loadV - loadY) - routeV->penalty;
    }

    if (costSuppU + costSuppV > -MY_EPSILON)
        return false;
    if (nodeX->isDepot || nodeY->isDepot || nodeY == nodeU->prev || nodeU == nodeY || nodeX == nodeV || nodeV == nodeX->next)
        return false;

    swapNode(nodeU, nodeV);
    swapNode(nodeX, nodeY);
    nbMoves++; // Increment move counter before updating route data
    searchCompleted = false;
    updateRouteData(routeU);
    if (!intraRouteMove)
        updateRouteData(routeV);
    return true;
}

bool LocalSearch::move7()
{
    if (nodeU->position > nodeV->position)
        return false;

    double cost = params.timeCost[nodeUIndex][nodeVIndex] + params.timeCost[nodeXIndex][nodeYIndex] - params.timeCost[nodeUIndex][nodeXIndex] - params.timeCost[nodeVIndex][nodeYIndex] + nodeV->cumulatedReversalDistance - nodeX->cumulatedReversalDistance;

    if (cost > -MY_EPSILON)
        return false;
    if (nodeU->next == nodeV)
        return false;

    Node *nodeNum = nodeX->next;
    nodeX->prev = nodeNum;
    nodeX->next = nodeY;

    while (nodeNum != nodeV)
    {
        Node *temp = nodeNum->next;
        nodeNum->next = nodeNum->prev;
        nodeNum->prev = temp;
        nodeNum = temp;
    }

    nodeV->next = nodeV->prev;
    nodeV->prev = nodeU;
    nodeU->next = nodeV;
    nodeY->prev = nodeX;

    nbMoves++; // Increment move counter before updating route data
    searchCompleted = false;
    updateRouteData(routeU);
    return true;
}

bool LocalSearch::move8()
{
    double cost = params.timeCost[nodeUIndex][nodeVIndex] + params.timeCost[nodeXIndex][nodeYIndex] - params.timeCost[nodeUIndex][nodeXIndex] - params.timeCost[nodeVIndex][nodeYIndex] + nodeV->cumulatedReversalDistance + routeU->reversalDistance - nodeX->cumulatedReversalDistance - routeU->penalty - routeV->penalty;

    // Early move pruning to save CPU time. Guarantees that this move cannot improve without checking additional (load, duration...) constraints
    if (cost >= 0)
        return false;

    cost += penaltyExcessDuration(nodeU->cumulatedTime + nodeV->cumulatedTime + nodeV->cumulatedReversalDistance + params.timeCost[nodeUIndex][nodeVIndex]) + penaltyExcessDuration(routeU->duration - nodeU->cumulatedTime - params.timeCost[nodeUIndex][nodeXIndex] + routeU->reversalDistance - nodeX->cumulatedReversalDistance + routeV->duration - nodeV->cumulatedTime - params.timeCost[nodeVIndex][nodeYIndex] + params.timeCost[nodeXIndex][nodeYIndex]) + penaltyExcessLoad(nodeU->cumulatedLoad + nodeV->cumulatedLoad) + penaltyExcessLoad(routeU->load + routeV->load - nodeU->cumulatedLoad - nodeV->cumulatedLoad);

    if (cost > -MY_EPSILON)
        return false;

    Node *depotU = routeU->depot;
    Node *depotV = routeV->depot;
    Node *depotUFin = routeU->depot->prev;
    Node *depotVFin = routeV->depot->prev;
    Node *depotVSuiv = depotV->next;

    Node *temp;
    Node *xx = nodeX;
    Node *vv = nodeV;

    while (!xx->isDepot)
    {
        temp = xx->next;
        xx->next = xx->prev;
        xx->prev = temp;
        xx->route = routeV;
        xx = temp;
    }

    while (!vv->isDepot)
    {
        temp = vv->prev;
        vv->prev = vv->next;
        vv->next = temp;
        vv->route = routeU;
        vv = temp;
    }

    nodeU->next = nodeV;
    nodeV->prev = nodeU;
    nodeX->next = nodeY;
    nodeY->prev = nodeX;

    if (nodeX->isDepot)
    {
        depotUFin->next = depotU;
        depotUFin->prev = depotVSuiv;
        depotUFin->prev->next = depotUFin;
        depotV->next = nodeY;
        nodeY->prev = depotV;
    }
    else if (nodeV->isDepot)
    {
        depotV->next = depotUFin->prev;
        depotV->next->prev = depotV;
        depotV->prev = depotVFin;
        depotUFin->prev = nodeU;
        nodeU->next = depotUFin;
    }
    else
    {
        depotV->next = depotUFin->prev;
        depotV->next->prev = depotV;
        depotUFin->prev = depotVSuiv;
        depotUFin->prev->next = depotUFin;
    }

    nbMoves++; // Increment move counter before updating route data
    searchCompleted = false;
    updateRouteData(routeU);
    updateRouteData(routeV);
    return true;
}

bool LocalSearch::move9()
{
    double cost = params.timeCost[nodeUIndex][nodeYIndex] + params.timeCost[nodeVIndex][nodeXIndex] - params.timeCost[nodeUIndex][nodeXIndex] - params.timeCost[nodeVIndex][nodeYIndex] - routeU->penalty - routeV->penalty;

    // Early move pruning to save CPU time. Guarantees that this move cannot improve without checking additional (load, duration...) constraints
    if (cost >= 0)
        return false;

    cost += penaltyExcessDuration(nodeU->cumulatedTime + routeV->duration - nodeV->cumulatedTime - params.timeCost[nodeVIndex][nodeYIndex] + params.timeCost[nodeUIndex][nodeYIndex]) + penaltyExcessDuration(routeU->duration - nodeU->cumulatedTime - params.timeCost[nodeUIndex][nodeXIndex] + nodeV->cumulatedTime + params.timeCost[nodeVIndex][nodeXIndex]) + penaltyExcessLoad(nodeU->cumulatedLoad + routeV->load - nodeV->cumulatedLoad) + penaltyExcessLoad(nodeV->cumulatedLoad + routeU->load - nodeU->cumulatedLoad);

    if (cost > -MY_EPSILON)
        return false;

    Node *depotU = routeU->depot;
    Node *depotV = routeV->depot;
    Node *depotUFin = depotU->prev;
    Node *depotVFin = depotV->prev;
    Node *depotUpred = depotUFin->prev;

    Node *count = nodeY;
    while (!count->isDepot)
    {
        count->route = routeU;
        count = count->next;
    }

    count = nodeX;
    while (!count->isDepot)
    {
        count->route = routeV;
        count = count->next;
    }

    nodeU->next = nodeY;
    nodeY->prev = nodeU;
    nodeV->next = nodeX;
    nodeX->prev = nodeV;

    if (nodeX->isDepot)
    {
        depotUFin->prev = depotVFin->prev;
        depotUFin->prev->next = depotUFin;
        nodeV->next = depotVFin;
        depotVFin->prev = nodeV;
    }
    else
    {
        depotUFin->prev = depotVFin->prev;
        depotUFin->prev->next = depotUFin;
        depotVFin->prev = depotUpred;
        depotVFin->prev->next = depotVFin;
    }

    nbMoves++; // Increment move counter before updating route data
    searchCompleted = false;
    updateRouteData(routeU);
    updateRouteData(routeV);
    return true;
}

bool LocalSearch::swapStar()
{
    SwapStarElement myBestSwapStar;

    // Preprocessing insertion costs
    preprocessInsertions(routeU, routeV);
    preprocessInsertions(routeV, routeU);

    // Evaluating the moves
    for (nodeU = routeU->depot->next; !nodeU->isDepot; nodeU = nodeU->next)
    {
        for (nodeV = routeV->depot->next; !nodeV->isDepot; nodeV = nodeV->next)
        {
            double deltaPenRouteU = penaltyExcessLoad(routeU->load + params.cli[nodeV->cour].demand - params.cli[nodeU->cour].demand) - routeU->penalty;
            double deltaPenRouteV = penaltyExcessLoad(routeV->load + params.cli[nodeU->cour].demand - params.cli[nodeV->cour].demand) - routeV->penalty;

            // Quick filter: possibly early elimination of many SWAP* due to the capacity constraints/penalties and bounds on insertion costs
            if (deltaPenRouteU + nodeU->deltaRemoval + deltaPenRouteV + nodeV->deltaRemoval <= 0)
            {
                SwapStarElement mySwapStar;
                mySwapStar.U = nodeU;
                mySwapStar.V = nodeV;

                // Evaluate best reinsertion cost of U in the route of V where V has been removed
                double extraV = getCheapestInsertSimultRemoval(nodeU, nodeV, mySwapStar.bestPositionU);

                // Evaluate best reinsertion cost of V in the route of U where U has been removed
                double extraU = getCheapestInsertSimultRemoval(nodeV, nodeU, mySwapStar.bestPositionV);

                // Evaluating final cost
                mySwapStar.moveCost = deltaPenRouteU + nodeU->deltaRemoval + extraU + deltaPenRouteV + nodeV->deltaRemoval + extraV + penaltyExcessDuration(routeU->duration + nodeU->deltaRemoval + extraU + params.cli[nodeV->cour].serviceDuration - params.cli[nodeU->cour].serviceDuration) + penaltyExcessDuration(routeV->duration + nodeV->deltaRemoval + extraV - params.cli[nodeV->cour].serviceDuration + params.cli[nodeU->cour].serviceDuration);

                if (mySwapStar.moveCost < myBestSwapStar.moveCost)
                    myBestSwapStar = mySwapStar;
            }
        }
    }

    // Including RELOCATE from nodeU towards routeV (costs nothing to include in the evaluation at this step since we already have the best insertion location)
    // Moreover, since the granularity criterion is different, this can lead to different improving moves
    for (nodeU = routeU->depot->next; !nodeU->isDepot; nodeU = nodeU->next)
    {
        SwapStarElement mySwapStar;
        mySwapStar.U = nodeU;
        mySwapStar.bestPositionU = bestInsertClient[routeV->cour][nodeU->cour].bestLocation[0];
        double deltaDistRouteU = params.timeCost[nodeU->prev->cour][nodeU->next->cour] - params.timeCost[nodeU->prev->cour][nodeU->cour] - params.timeCost[nodeU->cour][nodeU->next->cour];
        double deltaDistRouteV = bestInsertClient[routeV->cour][nodeU->cour].bestCost[0];
        mySwapStar.moveCost = deltaDistRouteU + deltaDistRouteV + penaltyExcessLoad(routeU->load - params.cli[nodeU->cour].demand) - routeU->penalty + penaltyExcessLoad(routeV->load + params.cli[nodeU->cour].demand) - routeV->penalty + penaltyExcessDuration(routeU->duration + deltaDistRouteU - params.cli[nodeU->cour].serviceDuration) + penaltyExcessDuration(routeV->duration + deltaDistRouteV + params.cli[nodeU->cour].serviceDuration);

        if (mySwapStar.moveCost < myBestSwapStar.moveCost)
            myBestSwapStar = mySwapStar;
    }

    // Including RELOCATE from nodeV towards routeU
    for (nodeV = routeV->depot->next; !nodeV->isDepot; nodeV = nodeV->next)
    {
        SwapStarElement mySwapStar;
        mySwapStar.V = nodeV;
        mySwapStar.bestPositionV = bestInsertClient[routeU->cour][nodeV->cour].bestLocation[0];
        double deltaDistRouteU = bestInsertClient[routeU->cour][nodeV->cour].bestCost[0];
        double deltaDistRouteV = params.timeCost[nodeV->prev->cour][nodeV->next->cour] - params.timeCost[nodeV->prev->cour][nodeV->cour] - params.timeCost[nodeV->cour][nodeV->next->cour];
        mySwapStar.moveCost = deltaDistRouteU + deltaDistRouteV + penaltyExcessLoad(routeU->load + params.cli[nodeV->cour].demand) - routeU->penalty + penaltyExcessLoad(routeV->load - params.cli[nodeV->cour].demand) - routeV->penalty + penaltyExcessDuration(routeU->duration + deltaDistRouteU + params.cli[nodeV->cour].serviceDuration) + penaltyExcessDuration(routeV->duration + deltaDistRouteV - params.cli[nodeV->cour].serviceDuration);

        if (mySwapStar.moveCost < myBestSwapStar.moveCost)
            myBestSwapStar = mySwapStar;
    }

    if (myBestSwapStar.moveCost > -MY_EPSILON)
        return false;

    // Applying the best move in case of improvement
    if (myBestSwapStar.bestPositionU != NULL)
        insertNode(myBestSwapStar.U, myBestSwapStar.bestPositionU);
    if (myBestSwapStar.bestPositionV != NULL)
        insertNode(myBestSwapStar.V, myBestSwapStar.bestPositionV);
    nbMoves++; // Increment move counter before updating route data
    searchCompleted = false;
    updateRouteData(routeU);
    updateRouteData(routeV);
    return true;
}

double LocalSearch::getCheapestInsertSimultRemoval(Node *U, Node *V, Node *&bestPosition)
{
    ThreeBestInsert *myBestInsert = &bestInsertClient[V->route->cour][U->cour];
    bool found = false;

    // Find best insertion in the route such that V is not next or pred (can only belong to the top three locations)
    bestPosition = myBestInsert->bestLocation[0];
    double bestCost = myBestInsert->bestCost[0];
    found = (bestPosition != V && bestPosition->next != V);
    if (!found && myBestInsert->bestLocation[1] != NULL)
    {
        bestPosition = myBestInsert->bestLocation[1];
        bestCost = myBestInsert->bestCost[1];
        found = (bestPosition != V && bestPosition->next != V);
        if (!found && myBestInsert->bestLocation[2] != NULL)
        {
            bestPosition = myBestInsert->bestLocation[2];
            bestCost = myBestInsert->bestCost[2];
            found = true;
        }
    }

    // Compute insertion in the place of V
    double deltaCost = params.timeCost[V->prev->cour][U->cour] + params.timeCost[U->cour][V->next->cour] - params.timeCost[V->prev->cour][V->next->cour];
    if (!found || deltaCost < bestCost)
    {
        bestPosition = V->prev;
        bestCost = deltaCost;
    }

    return bestCost;
}

void LocalSearch::preprocessInsertions(Route *R1, Route *R2)
{
    for (Node *U = R1->depot->next; !U->isDepot; U = U->next)
    {
        // Performs the preprocessing
        U->deltaRemoval = params.timeCost[U->prev->cour][U->next->cour] - params.timeCost[U->prev->cour][U->cour] - params.timeCost[U->cour][U->next->cour];
        if (R2->whenLastModified > bestInsertClient[R2->cour][U->cour].whenLastCalculated)
        {
            bestInsertClient[R2->cour][U->cour].reset();
            bestInsertClient[R2->cour][U->cour].whenLastCalculated = nbMoves;
            bestInsertClient[R2->cour][U->cour].bestCost[0] = params.timeCost[0][U->cour] + params.timeCost[U->cour][R2->depot->next->cour] - params.timeCost[0][R2->depot->next->cour];
            bestInsertClient[R2->cour][U->cour].bestLocation[0] = R2->depot;
            for (Node *V = R2->depot->next; !V->isDepot; V = V->next)
            {
                double deltaCost = params.timeCost[V->cour][U->cour] + params.timeCost[U->cour][V->next->cour] - params.timeCost[V->cour][V->next->cour];
                bestInsertClient[R2->cour][U->cour].compareAndAdd(deltaCost, V);
            }
        }
    }
}

void LocalSearch::insertNode(Node *U, Node *V)
{
    U->prev->next = U->next;
    U->next->prev = U->prev;
    V->next->prev = U;
    U->prev = V;
    U->next = V->next;
    V->next = U;
    U->route = V->route;
}

void LocalSearch::swapNode(Node *U, Node *V)
{
    Node *myVPred = V->prev;
    Node *myVSuiv = V->next;
    Node *myUPred = U->prev;
    Node *myUSuiv = U->next;
    Route *myRouteU = U->route;
    Route *myRouteV = V->route;

    myUPred->next = V;
    myUSuiv->prev = V;
    myVPred->next = U;
    myVSuiv->prev = U;

    U->prev = myVPred;
    U->next = myVSuiv;
    V->prev = myUPred;
    V->next = myUSuiv;

    U->route = myRouteV;
    V->route = myRouteU;
}

void LocalSearch::updateRouteData(Route *myRoute)
{
    int myplace = 0;
    double myload = 0.;
    double mytime = 0.;
    double myReversalDistance = 0.;
    double cumulatedX = 0.;
    double cumulatedY = 0.;

    Node *mynode = myRoute->depot;
    mynode->position = 0;
    mynode->cumulatedLoad = 0.;
    mynode->cumulatedTime = 0.;
    mynode->cumulatedReversalDistance = 0.;

    bool firstIt = true;
    while (!mynode->isDepot || firstIt)
    {
        mynode = mynode->next;
        myplace++;
        mynode->position = myplace;
        myload += params.cli[mynode->cour].demand;
        mytime += params.timeCost[mynode->prev->cour][mynode->cour] + params.cli[mynode->cour].serviceDuration;
        myReversalDistance += params.timeCost[mynode->cour][mynode->prev->cour] - params.timeCost[mynode->prev->cour][mynode->cour];
        mynode->cumulatedLoad = myload;
        mynode->cumulatedTime = mytime;
        mynode->cumulatedReversalDistance = myReversalDistance;
        if (!mynode->isDepot)
        {
            cumulatedX += params.cli[mynode->cour].coordX;
            cumulatedY += params.cli[mynode->cour].coordY;
            if (firstIt)
                myRoute->sector.initialize(params.cli[mynode->cour].polarAngle);
            else
                myRoute->sector.extend(params.cli[mynode->cour].polarAngle);
        }
        firstIt = false;
    }

    myRoute->duration = mytime;
    myRoute->load = myload;
    myRoute->penalty = penaltyExcessDuration(mytime) + penaltyExcessLoad(myload);
    myRoute->nbCustomers = myplace - 1;
    myRoute->reversalDistance = myReversalDistance;
    // Remember "when" this route has been last modified (will be used to filter unnecessary move evaluations)
    myRoute->whenLastModified = nbMoves;

    if (myRoute->nbCustomers == 0)
    {
        myRoute->polarAngleBarycenter = 1.e30;
        emptyRoutes.insert(myRoute->cour);
    }
    else
    {
        myRoute->polarAngleBarycenter = atan2(cumulatedY / (double)myRoute->nbCustomers - params.cli[0].coordY, cumulatedX / (double)myRoute->nbCustomers - params.cli[0].coordX);
        emptyRoutes.erase(myRoute->cour);
    }
}

void LocalSearch::loadIndividual(const Individual &indiv)
{
    emptyRoutes.clear();
    nbMoves = 0;
    for (int r = 0; r < params.nbVehicles; r++)
    {
        Node *myDepot = &depots[r];
        Node *myDepotFin = &depotsEnd[r];
        Route *myRoute = &routes[r];
        myDepot->prev = myDepotFin;
        myDepotFin->next = myDepot;
        if (!indiv.chromR[r].empty())
        {
            Node *myClient = &clients[indiv.chromR[r][0]];
            myClient->route = myRoute;
            myClient->prev = myDepot;
            myDepot->next = myClient;
            for (int i = 1; i < (int)indiv.chromR[r].size(); i++)
            {
                Node *myClientPred = myClient;
                myClient = &clients[indiv.chromR[r][i]];
                myClient->prev = myClientPred;
                myClientPred->next = myClient;
                myClient->route = myRoute;
            }
            myClient->next = myDepotFin;
            myDepotFin->prev = myClient;
        }
        else
        {
            myDepot->next = myDepotFin;
            myDepotFin->prev = myDepot;
        }
        updateRouteData(&routes[r]);
        routes[r].whenLastTestedSWAPStar = -1;
        // bookmark
        for (int i = 1; i <= params.nbClients; i++) // Initializing memory structures
            bestInsertClient[r][i].whenLastCalculated = -1;
    }
    // bookmark
    for (int i = 1; i <= params.nbClients; i++) // Initializing memory structures
        clients[i].whenLastTestedRI = -1;
}

void LocalSearch::exportIndividual(Individual &indiv)
{
    std::vector<std::pair<double, int>> routePolarAngles;
    for (int r = 0; r < params.nbVehicles; r++)
        routePolarAngles.push_back(std::pair<double, int>(routes[r].polarAngleBarycenter, r));
    std::sort(routePolarAngles.begin(), routePolarAngles.end()); // empty routes have a polar angle of 1.e30, and therefore will always appear at the end

    int pos = 0;
    // bookmark
    for (int r = 0; r < params.nbVehicles; r++)
    {
        indiv.chromR[r].clear();
        Node *node = depots[routePolarAngles[r].second].next;
        while (!node->isDepot)
        {
            indiv.chromT[pos] = node->cour;
            indiv.chromR[r].push_back(node->cour);
            node = node->next;
            pos++;
        }
    }

    indiv.evaluateCompleteCost(params);
}

LocalSearch::LocalSearch(Params &params) : params(params)
{
    clients = std::vector<Node>(params.nbClients + 1);
    routes = std::vector<Route>(params.nbVehicles);
    depots = std::vector<Node>(params.nbVehicles);
    depotsEnd = std::vector<Node>(params.nbVehicles);
    bestInsertClient = std::vector<std::vector<ThreeBestInsert>>(params.nbVehicles, std::vector<ThreeBestInsert>(params.nbClients + 1));

    // bookmark
    for (int i = 0; i <= params.nbClients; i++)
    {
        clients[i].cour = i;
        clients[i].isDepot = false;
    }
    // bookmark
    for (int i = 0; i < params.nbVehicles; i++)
    {
        routes[i].cour = i;
        routes[i].depot = &depots[i];
        depots[i].cour = 0;
        depots[i].isDepot = true;
        depots[i].route = &routes[i];
        depotsEnd[i].cour = 0;
        depotsEnd[i].isDepot = true;
        depotsEnd[i].route = &routes[i];
    }
    for (int i = 1; i <= params.nbClients; i++)
        orderNodes.push_back(i);
    for (int r = 0; r < params.nbVehicles; r++)
        orderRoutes.push_back(r);
}

#include "Params.h"

// The universal constructor for both executable and shared library
// When the executable is run from the commandline,
// it will first generate an CVRPLIB instance from .vrp file, then supply necessary information.
Params::Params(
    const std::vector<double> &x_coords,
    const std::vector<double> &y_coords,
    const std::vector<std::vector<double>> &dist_mtx,
    const std::vector<double> &service_time,
    const std::vector<double> &demands,
    double vehicleCapacity,
    double durationLimit,
    int nbVeh,
    bool isDurationConstraint,
    bool verbose,
    const AlgorithmParameters &ap)
    : ap(ap), isDurationConstraint(isDurationConstraint), nbVehicles(nbVeh), durationLimit(durationLimit),
      vehicleCapacity(vehicleCapacity), timeCost(dist_mtx), verbose(verbose)
{
    // This marks the starting time of the algorithm
    startTime = clock();

    nbClients = (int)demands.size() - 1; // Need to substract the depot from the number of nodes
    totalDemand = 0.;
    maxDemand = 0.;

    // Initialize RNG
    ran.seed(ap.seed);

    // check if valid coordinates are provided
    areCoordinatesProvided = (demands.size() == x_coords.size()) && (demands.size() == y_coords.size());

    cli = std::vector<Client>(nbClients + 1);
    // bookmark
    for (int i = 0; i <= nbClients; i++)
    {
        // If useSwapStar==false, x_coords and y_coords may be empty.
        if (ap.useSwapStar == 1 && areCoordinatesProvided)
        {
            cli[i].coordX = x_coords[i];
            cli[i].coordY = y_coords[i];
            cli[i].polarAngle = CircleSector::positive_mod(
                32768. * atan2(cli[i].coordY - cli[0].coordY, cli[i].coordX - cli[0].coordX) / PI);
        }
        else
        {
            cli[i].coordX = 0.0;
            cli[i].coordY = 0.0;
            cli[i].polarAngle = 0.0;
        }

        cli[i].serviceDuration = service_time[i];
        cli[i].demand = demands[i];
        if (cli[i].demand > maxDemand)
            maxDemand = cli[i].demand;
        totalDemand += cli[i].demand;
    }

    if (verbose && ap.useSwapStar == 1 && !areCoordinatesProvided)
        std::cout << "----- NO COORDINATES HAVE BEEN PROVIDED, SWAP* NEIGHBORHOOD WILL BE DEACTIVATED BY DEFAULT" << std::endl;

    // Default initialization if the number of vehicles has not been provided by the user
    if (nbVehicles == INT_MAX)
    {
        nbVehicles = (int)std::ceil(1.3 * totalDemand / vehicleCapacity) + 3; // Safety margin: 30% + 3 more vehicles than the trivial bin packing LB
        if (verbose)
            std::cout << "----- FLEET SIZE WAS NOT SPECIFIED: DEFAULT INITIALIZATION TO " << nbVehicles << " VEHICLES" << std::endl;
    }
    else
    {
        if (verbose)
            std::cout << "----- FLEET SIZE SPECIFIED: SET TO " << nbVehicles << " VEHICLES" << std::endl;
    }

    // Calculation of the maximum distance
    maxDist = 0.;
    // bookmark
    for (int i = 0; i <= nbClients; i++)
        for (int j = 0; j <= nbClients; j++)
            if (timeCost[i][j] > maxDist)
                maxDist = timeCost[i][j];

    // Calculation of the correlated vertices for each customer (for the granular restriction)
    correlatedVertices = std::vector<std::vector<int>>(nbClients + 1);
    std::vector<std::set<int>> setCorrelatedVertices = std::vector<std::set<int>>(nbClients + 1);
    std::vector<std::pair<double, int>> orderProximity;
    for (int i = 1; i <= nbClients; i++)
    {
        orderProximity.clear();
        for (int j = 1; j <= nbClients; j++)
            if (i != j)
                orderProximity.emplace_back(timeCost[i][j], j);
        std::sort(orderProximity.begin(), orderProximity.end());

        for (int j = 0; j < std::min<int>(ap.nbGranular, nbClients - 1); j++)
        {
            // If i is correlated with j, then j should be correlated with i
            setCorrelatedVertices[i].insert(orderProximity[j].second);
            setCorrelatedVertices[orderProximity[j].second].insert(i);
        }
    }

    // Filling the vector of correlated vertices
    // bookmark
    for (int i = 1; i <= nbClients; i++)
        for (int x : setCorrelatedVertices[i])
            correlatedVertices[i].push_back(x);

    // Safeguards to avoid possible numerical instability in case of instances containing arbitrarily small or large numerical values
    if (maxDist < 0.1 || maxDist > 100000)
        throw std::string(
            "The distances are of very small or large scale. This could impact numerical stability. Please rescale the dataset and run again.");
    if (maxDemand < 0.1 || maxDemand > 100000)
        throw std::string(
            "The demand quantities are of very small or large scale. This could impact numerical stability. Please rescale the dataset and run again.");
    if (nbVehicles < std::ceil(totalDemand / vehicleCapacity))
        throw std::string("Fleet size is insufficient to service the considered clients.");

    // A reasonable scale for the initial values of the penalties
    penaltyDuration = 1;
    penaltyCapacity = std::max<double>(0.1, std::min<double>(1000., maxDist / maxDemand));

    if (verbose)
        std::cout << "----- INSTANCE SUCCESSFULLY LOADED WITH " << nbClients << " CLIENTS AND " << nbVehicles << " VEHICLES" << std::endl;
}

#include "Population.h"

void Population::generatePopulation()
{
    if (params.verbose)
        std::cout << "----- BUILDING INITIAL POPULATION" << std::endl;
    for (int i = 0; i < 4 * params.ap.mu && (i == 0 || params.ap.timeLimit == 0 || (double)(clock() - params.startTime) < params.ap.timeLimit); i++)
    {
        Individual randomIndiv(params);
        split.generalSplit(randomIndiv, params.nbVehicles);
        localSearch.run(randomIndiv, params.penaltyCapacity, params.penaltyDuration);
        addIndividual(randomIndiv, true);
        if (!randomIndiv.eval.isFeasible && params.ran() % 2 == 0) // Repair half of the solutions in case of infeasibility
        {
            localSearch.run(randomIndiv, params.penaltyCapacity * 10., params.penaltyDuration * 10.);
            if (randomIndiv.eval.isFeasible)
                addIndividual(randomIndiv, false);
        }
    }
}

bool Population::addIndividual(const Individual &indiv, bool updateFeasible)
{
    if (updateFeasible)
    {
        listFeasibilityLoad.push_back(indiv.eval.capacityExcess < MY_EPSILON);
        listFeasibilityDuration.push_back(indiv.eval.durationExcess < MY_EPSILON);
        listFeasibilityLoad.pop_front();
        listFeasibilityDuration.pop_front();
    }

    // Find the adequate subpopulation in relation to the individual feasibility
    SubPopulation &subpop = (indiv.eval.isFeasible) ? feasibleSubpop : infeasibleSubpop;

    // Create a copy of the individual and updade the proximity structures calculating inter-individual distances
    Individual *myIndividual = new Individual(indiv);
    for (Individual *myIndividual2 : subpop)
    {
        double myDistance = brokenPairsDistance(*myIndividual, *myIndividual2);
        myIndividual2->indivsPerProximity.insert({myDistance, myIndividual});
        myIndividual->indivsPerProximity.insert({myDistance, myIndividual2});
    }

    // Identify the correct location in the subpopulation and insert the individual
    int place = (int)subpop.size();
    while (place > 0 && subpop[place - 1]->eval.penalizedCost > indiv.eval.penalizedCost - MY_EPSILON)
        place--;
    subpop.emplace(subpop.begin() + place, myIndividual);

    // Trigger a survivor selection if the maximimum subpopulation size is exceeded
    if ((int)subpop.size() > params.ap.mu + params.ap.lambda)
        while ((int)subpop.size() > params.ap.mu)
            removeWorstBiasedFitness(subpop);

    // Track best solution
    if (indiv.eval.isFeasible && indiv.eval.penalizedCost < bestSolutionRestart.eval.penalizedCost - MY_EPSILON)
    {
        bestSolutionRestart = indiv; // Copy
        if (indiv.eval.penalizedCost < bestSolutionOverall.eval.penalizedCost - MY_EPSILON)
        {
            bestSolutionOverall = indiv;
            searchProgress.push_back({clock() - params.startTime, bestSolutionOverall.eval.penalizedCost});
        }
        return true;
    }
    else
        return false;
}

__global__ void updateBiasedFitnesses_kernel(int pop_size, int *ranking_second, int params_ap_nbElite, double *pop_ranking_second_biasedFitness)
{
    int tid = threadIdx.x;

    if (tid < (int)pop_size)
    {
        double divRank = (double)tid / (double)(pop_size - 1); // Ranking from 0 to 1
        double fitRank = (double)ranking_second[tid] / (double)(pop_size - 1);
        if ((int)pop_size <= params_ap_nbElite) // Elite individuals cannot be smaller than population size
            pop_ranking_second_biasedFitness[tid] = fitRank;
        else
            pop_ranking_second_biasedFitness[tid] = fitRank + (1.0 - (double)params_ap_nbElite / (double)pop_size) * divRank;
    }
}

void Population::updateBiasedFitnesses(SubPopulation &pop)
{
    // Ranking the individuals based on their diversity contribution (decreasing order of distance)
    std::vector<std::pair<double, int>> ranking;
    for (int i = 0; i < (int)pop.size(); i++)
        ranking.push_back({-averageBrokenPairsDistanceClosest(*pop[i], params.ap.nbClose), i});
    std::sort(ranking.begin(), ranking.end());

    // Updating the biased fitness values
    if (pop.size() == 1)
        pop[0]->biasedFitness = 0;
    else
    {
        // bookmarkimp
        int pop_size = (int)pop.size();
        int params_ap_nbElite = params.ap.nbElite;
        int *ranking_second, *ranking_second2=(int *)malloc(sizeof(int) * pop_size);
        double *pop_ranking_second_biasedFitness, *pop_ranking_second_biasedFitness2=(double *)malloc(sizeof(double) * pop_size);

        for (int i = 0; i < pop_size; i++)
        {
            ranking_second2[i] = ranking[i].second;
        }

        cudaMalloc((void **)&ranking_second, pop_size * sizeof(int));
        cudaMalloc((void **)&pop_ranking_second_biasedFitness, pop_size * sizeof(int));

        cudaMemcpy(ranking_second, ranking_second2, pop_size * sizeof(int), cudaMemcpyHostToDevice);

        updateBiasedFitnesses_kernel<<<BLOCKS, NUM_THREADS>>>(pop_size, ranking_second, params_ap_nbElite, pop_ranking_second_biasedFitness);

        cudaMemcpy(pop_ranking_second_biasedFitness2, pop_ranking_second_biasedFitness, pop_size * sizeof(int), cudaMemcpyDeviceToHost);

        for (int i = 0; i < pop_size; i++)
        {
            pop[ranking[i].second]->biasedFitness = pop_ranking_second_biasedFitness2[i];
        }

        /*
        for (int i = 0; i < (int)pop.size(); i++)
        {
            double divRank = (double)i / (double)(pop.size() - 1); // Ranking from 0 to 1
            double fitRank = (double)ranking[i].second / (double)(pop.size() - 1);
            if ((int)pop.size() <= params.ap.nbElite) // Elite individuals cannot be smaller than population size
                pop[ranking[i].second]->biasedFitness = fitRank;
            else
                pop[ranking[i].second]->biasedFitness = fitRank + (1.0 - (double)params.ap.nbElite / (double)pop.size()) * divRank;
        }*/
    }
}

void Population::removeWorstBiasedFitness(SubPopulation &pop)
{
    updateBiasedFitnesses(pop);
    if (pop.size() <= 1)
        throw std::string("Eliminating the best individual: this should not occur in HGS");

    Individual *worstIndividual = NULL;
    int worstIndividualPosition = -1;
    bool isWorstIndividualClone = false;
    double worstIndividualBiasedFitness = -1.e30;
    for (int i = 1; i < (int)pop.size(); i++)
    {
        bool isClone = (averageBrokenPairsDistanceClosest(*pop[i], 1) < MY_EPSILON); // A distance equal to 0 indicates that a clone exists
        if ((isClone && !isWorstIndividualClone) || (isClone == isWorstIndividualClone && pop[i]->biasedFitness > worstIndividualBiasedFitness))
        {
            worstIndividualBiasedFitness = pop[i]->biasedFitness;
            isWorstIndividualClone = isClone;
            worstIndividualPosition = i;
            worstIndividual = pop[i];
        }
    }

    // Removing the individual from the population and freeing memory
    pop.erase(pop.begin() + worstIndividualPosition);

    // Cleaning its distances from the other individuals in the population
    for (Individual *indiv2 : pop)
    {
        auto it = indiv2->indivsPerProximity.begin();
        while (it->second != worstIndividual)
            ++it;
        indiv2->indivsPerProximity.erase(it);
    }

    // Freeing memory
    delete worstIndividual;
}

void Population::restart()
{
    if (params.verbose)
        std::cout << "----- RESET: CREATING A NEW POPULATION -----" << std::endl;
    for (Individual *indiv : feasibleSubpop)
        delete indiv;
    for (Individual *indiv : infeasibleSubpop)
        delete indiv;
    feasibleSubpop.clear();
    infeasibleSubpop.clear();
    bestSolutionRestart = Individual(params);
    generatePopulation();
}

void Population::managePenalties()
{
    // Setting some bounds [0.1,100000] to the penalty values for safety
    double fractionFeasibleLoad = (double)std::count(listFeasibilityLoad.begin(), listFeasibilityLoad.end(), true) / (double)listFeasibilityLoad.size();
    if (fractionFeasibleLoad < params.ap.targetFeasible - 0.05 && params.penaltyCapacity < 100000.)
        params.penaltyCapacity = std::min<double>(params.penaltyCapacity * params.ap.penaltyIncrease, 100000.);
    else if (fractionFeasibleLoad > params.ap.targetFeasible + 0.05 && params.penaltyCapacity > 0.1)
        params.penaltyCapacity = std::max<double>(params.penaltyCapacity * params.ap.penaltyDecrease, 0.1);

    // Setting some bounds [0.1,100000] to the penalty values for safety
    double fractionFeasibleDuration = (double)std::count(listFeasibilityDuration.begin(), listFeasibilityDuration.end(), true) / (double)listFeasibilityDuration.size();
    if (fractionFeasibleDuration < params.ap.targetFeasible - 0.05 && params.penaltyDuration < 100000.)
        params.penaltyDuration = std::min<double>(params.penaltyDuration * params.ap.penaltyIncrease, 100000.);
    else if (fractionFeasibleDuration > params.ap.targetFeasible + 0.05 && params.penaltyDuration > 0.1)
        params.penaltyDuration = std::max<double>(params.penaltyDuration * params.ap.penaltyDecrease, 0.1);

    // Update the evaluations
    // bookmarkimp
    for (int i = 0; i < (int)infeasibleSubpop.size(); i++)
        infeasibleSubpop[i]->eval.penalizedCost = infeasibleSubpop[i]->eval.distance + params.penaltyCapacity * infeasibleSubpop[i]->eval.capacityExcess + params.penaltyDuration * infeasibleSubpop[i]->eval.durationExcess;

    // If needed, reorder the individuals in the infeasible subpopulation since the penalty values have changed (simple bubble sort for the sake of simplicity)
    // bookmarkimp
    for (int i = 0; i < (int)infeasibleSubpop.size(); i++)
    {
        for (int j = 0; j < (int)infeasibleSubpop.size() - i - 1; j++)
        {
            if (infeasibleSubpop[j]->eval.penalizedCost > infeasibleSubpop[j + 1]->eval.penalizedCost + MY_EPSILON)
            {
                Individual *indiv = infeasibleSubpop[j];
                infeasibleSubpop[j] = infeasibleSubpop[j + 1];
                infeasibleSubpop[j + 1] = indiv;
            }
        }
    }
}

const Individual &Population::getBinaryTournament()
{
    // Picking two individuals with uniform distribution over the union of the feasible and infeasible subpopulations
    std::uniform_int_distribution<> distr(0, feasibleSubpop.size() + infeasibleSubpop.size() - 1);
    int place1 = distr(params.ran);
    int place2 = distr(params.ran);
    Individual *indiv1 = (place1 >= (int)feasibleSubpop.size()) ? infeasibleSubpop[place1 - feasibleSubpop.size()] : feasibleSubpop[place1];
    Individual *indiv2 = (place2 >= (int)feasibleSubpop.size()) ? infeasibleSubpop[place2 - feasibleSubpop.size()] : feasibleSubpop[place2];

    // Keeping the best of the two in terms of biased fitness
    updateBiasedFitnesses(feasibleSubpop);
    updateBiasedFitnesses(infeasibleSubpop);
    if (indiv1->biasedFitness < indiv2->biasedFitness)
        return *indiv1;
    else
        return *indiv2;
}

const Individual *Population::getBestFeasible()
{
    if (!feasibleSubpop.empty())
        return feasibleSubpop[0];
    else
        return NULL;
}

const Individual *Population::getBestInfeasible()
{
    if (!infeasibleSubpop.empty())
        return infeasibleSubpop[0];
    else
        return NULL;
}

const Individual *Population::getBestFound()
{
    if (bestSolutionOverall.eval.penalizedCost < 1.e29)
        return &bestSolutionOverall;
    else
        return NULL;
}

void Population::printState(int nbIter, int nbIterNoImprovement)
{
    if (params.verbose)
    {
        std::printf("It %6d %6d | T(s) %.2f", nbIter, nbIterNoImprovement, (double)(clock() - params.startTime) / (double)CLOCKS_PER_SEC);

        if (getBestFeasible() != NULL)
            std::printf(" | Feas %zu %.2f %.2f", feasibleSubpop.size(), getBestFeasible()->eval.penalizedCost, getAverageCost(feasibleSubpop));
        else
            std::printf(" | NO-FEASIBLE");

        if (getBestInfeasible() != NULL)
            std::printf(" | Inf %zu %.2f %.2f", infeasibleSubpop.size(), getBestInfeasible()->eval.penalizedCost, getAverageCost(infeasibleSubpop));
        else
            std::printf(" | NO-INFEASIBLE");

        std::printf(" | Div %.2f %.2f", getDiversity(feasibleSubpop), getDiversity(infeasibleSubpop));
        std::printf(" | Feas %.2f %.2f", (double)std::count(listFeasibilityLoad.begin(), listFeasibilityLoad.end(), true) / (double)listFeasibilityLoad.size(), (double)std::count(listFeasibilityDuration.begin(), listFeasibilityDuration.end(), true) / (double)listFeasibilityDuration.size());
        std::printf(" | Pen %.2f %.2f", params.penaltyCapacity, params.penaltyDuration);
        std::cout << std::endl;
    }
}

double Population::brokenPairsDistance(const Individual &indiv1, const Individual &indiv2)
{
    int differences = 0;
    // bookmarkreduction
    for (int j = 1; j <= params.nbClients; j++)
    {
        if (indiv1.successors[j] != indiv2.successors[j] && indiv1.successors[j] != indiv2.predecessors[j])
            differences++;
        if (indiv1.predecessors[j] == 0 && indiv2.predecessors[j] != 0 && indiv2.successors[j] != 0)
            differences++;
    }
    return (double)differences / (double)params.nbClients;
}

double Population::averageBrokenPairsDistanceClosest(const Individual &indiv, int nbClosest)
{
    double result = 0.;
    int maxSize = std::min<int>(nbClosest, indiv.indivsPerProximity.size());
    auto it = indiv.indivsPerProximity.begin();
    // bookmark
    for (int i = 0; i < maxSize; i++)
    {
        result += it->first;
        ++it;
    }
    return result / (double)maxSize;
}

double Population::getDiversity(const SubPopulation &pop)
{
    double average = 0.;
    int size = std::min<int>(params.ap.mu, pop.size()); // Only monitoring the "mu" better solutions to avoid too much noise in the measurements
    // bookmark
    for (int i = 0; i < size; i++)
        average += averageBrokenPairsDistanceClosest(*pop[i], size);
    if (size > 0)
        return average / (double)size;
    else
        return -1.0;
}

double Population::getAverageCost(const SubPopulation &pop)
{
    double average = 0.;
    int size = std::min<int>(params.ap.mu, pop.size()); // Only monitoring the "mu" better solutions to avoid too much noise in the measurements
    // bookmark
    for (int i = 0; i < size; i++)
        average += pop[i]->eval.penalizedCost;
    if (size > 0)
        return average / (double)size;
    else
        return -1.0;
}

void Population::exportSearchProgress(std::string fileName, std::string instanceName)
{
    std::ofstream myfile(fileName);
    for (std::pair<clock_t, double> state : searchProgress)
        myfile << instanceName << ";" << params.ap.seed << ";" << state.second << ";" << (double)state.first / (double)CLOCKS_PER_SEC << std::endl;
}

void Population::exportCVRPLibFormat(const Individual &indiv, std::string fileName)
{
    std::ofstream myfile(fileName);
    if (myfile.is_open())
    {
        for (int k = 0; k < (int)indiv.chromR.size(); k++)
        {
            if (!indiv.chromR[k].empty())
            {
                myfile << "Route #" << k + 1 << ":"; // Route IDs start at 1 in the file format
                for (int i : indiv.chromR[k])
                    myfile << " " << i;
                myfile << std::endl;
            }
        }
        myfile << "Cost " << indiv.eval.penalizedCost << std::endl;
    }
    else
        std::cout << "----- IMPOSSIBLE TO OPEN: " << fileName << std::endl;
}

Population::Population(Params &params, Split &split, LocalSearch &localSearch) : params(params), split(split), localSearch(localSearch), bestSolutionRestart(params), bestSolutionOverall(params)
{
    listFeasibilityLoad = std::list<bool>(params.ap.nbIterPenaltyManagement, true);
    listFeasibilityDuration = std::list<bool>(params.ap.nbIterPenaltyManagement, true);
}

Population::~Population()
{
    // bookmark
    for (int i = 0; i < (int)feasibleSubpop.size(); i++)
        delete feasibleSubpop[i];
    for (int i = 0; i < (int)infeasibleSubpop.size(); i++)
        delete infeasibleSubpop[i];
}

#include "Split.h"

void Split::generalSplit(Individual &indiv, int nbMaxVehicles)
{
    // Do not apply Split with fewer vehicles than the trivial (LP) bin packing bound
    maxVehicles = std::max<int>(nbMaxVehicles, std::ceil(params.totalDemand / params.vehicleCapacity));

    // Initialization of the data structures for the linear split algorithms
    // Direct application of the code located at https://github.com/vidalt/Split-Library
    // bookmarkimp
    for (int i = 1; i <= params.nbClients; i++)
    {
        cliSplit[i].demand = params.cli[indiv.chromT[i - 1]].demand;
        cliSplit[i].serviceTime = params.cli[indiv.chromT[i - 1]].serviceDuration;
        cliSplit[i].d0_x = params.timeCost[0][indiv.chromT[i - 1]];
        cliSplit[i].dx_0 = params.timeCost[indiv.chromT[i - 1]][0];
        if (i < params.nbClients)
            cliSplit[i].dnext = params.timeCost[indiv.chromT[i - 1]][indiv.chromT[i]];
        else
            cliSplit[i].dnext = -1.e30;
        sumLoad[i] = sumLoad[i - 1] + cliSplit[i].demand;
        sumService[i] = sumService[i - 1] + cliSplit[i].serviceTime;
        sumDistance[i] = sumDistance[i - 1] + cliSplit[i - 1].dnext;
    }

    // We first try the simple split, and then the Split with limited fleet if this is not successful
    if (splitSimple(indiv) == 0)
        splitLF(indiv);

    // Build up the rest of the Individual structure
    indiv.evaluateCompleteCost(params);
}

int Split::splitSimple(Individual &indiv)
{
    // Reinitialize the potential structures
    potential[0][0] = 0;
    // bookmarkpartiallyimp
    for (int i = 1; i <= params.nbClients; i++)
        potential[0][i] = 1.e30;

    // MAIN ALGORITHM -- Simple Split using Bellman's algorithm in topological order
    // This code has been maintained as it is very simple and can be easily adapted to a variety of constraints, whereas the O(n) Split has a more restricted application scope
    if (params.isDurationConstraint)
    {
        for (int i = 0; i < params.nbClients; i++)
        {
            double load = 0.;
            double distance = 0.;
            double serviceDuration = 0.;
            // bookmarkreduction
            for (int j = i + 1; j <= params.nbClients && load <= 1.5 * params.vehicleCapacity; j++)
            {
                load += cliSplit[j].demand;
                serviceDuration += cliSplit[j].serviceTime;
                if (j == i + 1)
                    distance += cliSplit[j].d0_x;
                else
                    distance += cliSplit[j - 1].dnext;
                double cost = distance + cliSplit[j].dx_0 + params.penaltyCapacity * std::max<double>(load - params.vehicleCapacity, 0.) + params.penaltyDuration * std::max<double>(distance + cliSplit[j].dx_0 + serviceDuration - params.durationLimit, 0.);
                if (potential[0][i] + cost < potential[0][j])
                {
                    potential[0][j] = potential[0][i] + cost;
                    pred[0][j] = i;
                }
            }
        }
    }
    else
    {
        Trivial_Deque queue = Trivial_Deque(params.nbClients + 1, 0);
        for (int i = 1; i <= params.nbClients; i++)
        {
            // The front is the best predecessor for i
            potential[0][i] = propagate(queue.get_front(), i, 0);
            pred[0][i] = queue.get_front();

            if (i < params.nbClients)
            {
                // If i is not dominated by the last of the pile
                if (!dominates(queue.get_back(), i, 0))
                {
                    // then i will be inserted, need to remove whoever is dominated by i.
                    while (queue.size() > 0 && dominatesRight(queue.get_back(), i, 0))
                        queue.pop_back();
                    queue.push_back(i);
                }
                // Check iteratively if front is dominated by the next front
                while (queue.size() > 1 && propagate(queue.get_front(), i + 1, 0) > propagate(queue.get_next_front(), i + 1, 0) - MY_EPSILON)
                    queue.pop_front();
            }
        }
    }

    if (potential[0][params.nbClients] > 1.e29)
        throw std::string("ERROR : no Split solution has been propagated until the last node");

    // Filling the chromR structure
    // bookmark
    for (int k = params.nbVehicles - 1; k >= maxVehicles; k--)
        indiv.chromR[k].clear();

    int end = params.nbClients;
    for (int k = maxVehicles - 1; k >= 0; k--)
    {
        indiv.chromR[k].clear();
        int begin = pred[0][end];
        for (int ii = begin; ii < end; ii++)
            indiv.chromR[k].push_back(indiv.chromT[ii]);
        end = begin;
    }

    // Return OK in case the Split algorithm reached the beginning of the routes
    return (end == 0);
}

// Split for problems with limited fleet
int Split::splitLF(Individual &indiv)
{
    // Initialize the potential structures
    potential[0][0] = 0;
    // bookmark
    for (int k = 0; k <= maxVehicles; k++)
        for (int i = 1; i <= params.nbClients; i++)
            potential[k][i] = 1.e30;

    // MAIN ALGORITHM -- Simple Split using Bellman's algorithm in topological order
    // This code has been maintained as it is very simple and can be easily adapted to a variety of constraints, whereas the O(n) Split has a more restricted application scope
    if (params.isDurationConstraint)
    {
        for (int k = 0; k < maxVehicles; k++)
        {
            for (int i = k; i < params.nbClients && potential[k][i] < 1.e29; i++)
            {
                double load = 0.;
                double serviceDuration = 0.;
                double distance = 0.;
                for (int j = i + 1; j <= params.nbClients && load <= 1.5 * params.vehicleCapacity; j++) // Setting a maximum limit on load infeasibility to accelerate the algorithm
                {
                    load += cliSplit[j].demand;
                    serviceDuration += cliSplit[j].serviceTime;
                    if (j == i + 1)
                        distance += cliSplit[j].d0_x;
                    else
                        distance += cliSplit[j - 1].dnext;
                    double cost = distance + cliSplit[j].dx_0 + params.penaltyCapacity * std::max<double>(load - params.vehicleCapacity, 0.) + params.penaltyDuration * std::max<double>(distance + cliSplit[j].dx_0 + serviceDuration - params.durationLimit, 0.);
                    if (potential[k][i] + cost < potential[k + 1][j])
                    {
                        potential[k + 1][j] = potential[k][i] + cost;
                        pred[k + 1][j] = i;
                    }
                }
            }
        }
    }
    else // MAIN ALGORITHM -- Without duration constraints in O(n), from "Vidal, T. (2016). Split algorithm in O(n) for the capacitated vehicle routing problem. C&OR"
    {
        Trivial_Deque queue = Trivial_Deque(params.nbClients + 1, 0);
        for (int k = 0; k < maxVehicles; k++)
        {
            // in the Split problem there is always one feasible solution with k routes that reaches the index k in the tour.
            queue.reset(k);

            // The range of potentials < 1.29 is always an interval.
            // The size of the queue will stay >= 1 until we reach the end of this interval.
            for (int i = k + 1; i <= params.nbClients && queue.size() > 0; i++)
            {
                // The front is the best predecessor for i
                potential[k + 1][i] = propagate(queue.get_front(), i, k);
                pred[k + 1][i] = queue.get_front();

                if (i < params.nbClients)
                {
                    // If i is not dominated by the last of the pile
                    if (!dominates(queue.get_back(), i, k))
                    {
                        // then i will be inserted, need to remove whoever he dominates
                        while (queue.size() > 0 && dominatesRight(queue.get_back(), i, k))
                            queue.pop_back();
                        queue.push_back(i);
                    }

                    // Check iteratively if front is dominated by the next front
                    while (queue.size() > 1 && propagate(queue.get_front(), i + 1, k) > propagate(queue.get_next_front(), i + 1, k) - MY_EPSILON)
                        queue.pop_front();
                }
            }
        }
    }

    if (potential[maxVehicles][params.nbClients] > 1.e29)
        throw std::string("ERROR : no Split solution has been propagated until the last node");

    // It could be cheaper to use a smaller number of vehicles
    double minCost = potential[maxVehicles][params.nbClients];
    int nbRoutes = maxVehicles;
    for (int k = 1; k < maxVehicles; k++)
        if (potential[k][params.nbClients] < minCost)
        {
            minCost = potential[k][params.nbClients];
            nbRoutes = k;
        }

    // Filling the chromR structure
    // bookmark
    for (int k = params.nbVehicles - 1; k >= nbRoutes; k--)
        indiv.chromR[k].clear();

    int end = params.nbClients;
    for (int k = nbRoutes - 1; k >= 0; k--)
    {
        indiv.chromR[k].clear();
        int begin = pred[k + 1][end];
        for (int ii = begin; ii < end; ii++)
            indiv.chromR[k].push_back(indiv.chromT[ii]);
        end = begin;
    }

    // Return OK in case the Split algorithm reached the beginning of the routes
    return (end == 0);
}

Split::Split(const Params &params) : params(params)
{
    // Structures of the linear Split
    cliSplit = std::vector<ClientSplit>(params.nbClients + 1);
    sumDistance = std::vector<double>(params.nbClients + 1, 0.);
    sumLoad = std::vector<double>(params.nbClients + 1, 0.);
    sumService = std::vector<double>(params.nbClients + 1, 0.);
    potential = std::vector<std::vector<double>>(params.nbVehicles + 1, std::vector<double>(params.nbClients + 1, 1.e30));
    pred = std::vector<std::vector<int>>(params.nbVehicles + 1, std::vector<int>(params.nbClients + 1, 0));
}

#include "Genetic.h"
#include "commandline.h"
#include "LocalSearch.h"
#include "Split.h"
#include "InstanceCVRPLIB.h"
#include "AlgorithmParameters.h"
using namespace std;

int main(int argc, char *argv[])
{
    try
    {
        // Reading the arguments of the program
        CommandLine commandline(argc, argv);

        // Print all algorithm parameter values
        if (commandline.verbose)
            print_algorithm_parameters(commandline.ap);

        // Reading the data file and initializing some data structures
        if (commandline.verbose)
            std::cout << "----- READING INSTANCE: " << commandline.pathInstance << std::endl;
        InstanceCVRPLIB cvrp(commandline.pathInstance, commandline.isRoundingInteger);

        Params params(cvrp.x_coords, cvrp.y_coords, cvrp.dist_mtx, cvrp.service_time, cvrp.demands,
                      cvrp.vehicleCapacity, cvrp.durationLimit, commandline.nbVeh, cvrp.isDurationConstraint, commandline.verbose, commandline.ap);

        // Running HGS
        Genetic solver(params);
        solver.run();

        // Exporting the best solution
        if (solver.population.getBestFound() != NULL)
        {
            if (params.verbose)
                std::cout << "----- WRITING BEST SOLUTION IN : " << commandline.pathSolution << std::endl;
            solver.population.exportCVRPLibFormat(*solver.population.getBestFound(), commandline.pathSolution);
            solver.population.exportSearchProgress(commandline.pathSolution + ".PG.csv", commandline.pathInstance);
        }
    }
    catch (const string &e)
    {
        std::cout << "EXCEPTION | " << e << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cout << "EXCEPTION | " << e.what() << std::endl;
    }
    return 0;
}
