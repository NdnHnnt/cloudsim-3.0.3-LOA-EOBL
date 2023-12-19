package org.cloudbus.cloudsim.examples;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.cloudbus.cloudsim.Cloudlet;
import org.cloudbus.cloudsim.Vm;

public class LionOptimizationAlgorithm {
    private int populationSize;
    private List<Cloudlet> cloudletList;
    private List<Vm> vmList;
    private double mutationRate;
    private List<Integer> vmAllocation;
    private List<Integer> cloudletAllocation;
    private double localBestFitness;
    private List<Integer> bestVmAllocation;
    private List<Integer> bestCloudletAllocation;
    private double globalBestFitness;

    public LionOptimizationAlgorithm(int populationSize, double mutationRate, List<Cloudlet> cloudletList, List<Vm> vmList) {
        this.populationSize = populationSize;
        this.mutationRate = mutationRate;
        this.cloudletList = cloudletList;
        this.vmList = vmList;
        this.localBestFitness = 0;
        this.globalBestFitness = 0;
    }

    public Population initPopulation(int chromosomeLength, int dataCenterIterator) {
        Population population = new Population(this.populationSize, chromosomeLength, dataCenterIterator);
        return population;
    }

    public void evalPopulation(Population population, int dataCenterIterator, int cloudletIteration) {
        for (Individual individual : population.getIndividuals()) {
            double individualFitness = calcFitness(individual, dataCenterIterator, cloudletIteration);
            individual.setFitness(individualFitness);
        }
    }

    public double calcFitness(Individual individual, int dataCenterIterator, int cloudletIteration) {
        double totalExecutionTime = 0;
        double totalCost = 0;
        int iterator = 0;
        dataCenterIterator = dataCenterIterator - 1;

        for (int i = 0 + dataCenterIterator * 9 + cloudletIteration * 54; i < 9 + dataCenterIterator * 9 + cloudletIteration * 54; i++) {
            int gene = individual.getGene(iterator);
            double mips = calculateMips(gene % 9);

            totalExecutionTime += cloudletList.get(i).getCloudletLength() / mips;
            totalCost += calculateCost(vmList.get(gene % 9), cloudletList.get(i));
            iterator++;
        }

        // Modify the fitness calculation to prioritize makespan optimization
        double makespanFitness = calculateMakespanFitness(totalExecutionTime);
        double costFitness = calculateCostFitness(totalCost);
        double fitness = makespanFitness + costFitness;

        individual.setFitness(fitness);
        return fitness;
    }

    public void mutatePopulation(Population population, int dataCenterIterator) {
        Individual[] individuals = population.getIndividuals();

        for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
            Individual individual = individuals[populationIndex];

            for (int geneIndex = 0; geneIndex < individual.getCloudletLength(); geneIndex++) {
                if (Math.random() < this.mutationRate) {
                    // Randomly select a new VM gene that is not equal to the current gene
                    int currentGene = individual.getGene(geneIndex);
                    int newGene = generateRandomNonEqualVmgene(currentGene, dataCenterIterator);
                    individual.setGene(geneIndex, newGene);
                }
            }
        }
    }

    public Population applyLionBehavior(Population population, int dataCenterIterator) {
        Individual[]individuals = population.getIndividuals();

        for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
            Individual individual = individuals[populationIndex];

            for (int geneIndex = 0; geneIndex < individual.getCloudletLength(); geneIndex++) {
                if (Math.random() < this.mutationRate) {
                    int currentGene = individual.getGene(geneIndex);
                    int newGene = applyLionBehaviorToGene(currentGene, dataCenterIterator);
                    individual.setGene(geneIndex, newGene);
                }
            }
        }

        return population;
    }

    public int generateRandomNonEqualVmgene(int currentGene, int dataCenterIterator) {
        int minGene = dataCenterIterator * 9;
        int maxGene = (dataCenterIterator + 1) * 9 - 1;
        int newGene;

        // Generate a random VM gene until it is not equal to the current gene
        do {
            newGene = new Random().nextInt(maxGene - minGene + 1) + minGene;
        } while (newGene == currentGene);

        return newGene;
    }

    public int applyLionBehaviorToGene(int gene, int dataCenterIterator) {
        int minGene = dataCenterIterator * 9;
        int maxGene = (dataCenterIterator + 1) * 9 - 1;

        if (gene < minGene || gene > maxGene) {
            return generateRandomNonEqualVmgene(gene, dataCenterIterator);
        }

        int squirrelJump = 1; // Adjust the squirrel jump value based on the problem domain

        // Apply squirrel movement by adding/subtracting squirrelJump to the current gene
        int newGene = gene + squirrelJump;
        if (newGene > maxGene) {
            newGene = maxGene;
        } else if (newGene < minGene) {
            newGene = minGene;
        }

        return newGene;
    }

    private double calculateMips(int vmIndex) {
        double mips = 0;
        if (vmIndex % 9 == 0 || vmIndex % 9 == 3 || vmIndex % 9 == 6) {
            mips = 400;
        } else if (vmIndex % 9 == 1 || vmIndex % 9 == 4 || vmIndex % 9 == 7) {
            mips = 500;
        } else if (vmIndex % 9 == 2 || vmIndex % 9 == 5 || vmIndex % 9 == 8) {
            mips = 600;
        }
        return mips;
    }

    private double calculateCost(Vm vm, Cloudlet cloudlet) {
        double costPerMips = vm.getCostPerMips();
        double cloudletLength = cloudlet.getCloudletLength();
        double mips = vm.getMips();
        double executionTime = cloudletLength / mips;
        return costPerMips * executionTime;
    }

    private double calculateMakespanFitness(double totalExecutionTime) {
        // The higher the makespan, the lower the fitness
        return 1.0 / totalExecutionTime;
    }

    private double calculateCostFitness(double totalCost) {
        // The lower the cost, the higher the fitness
        return 1.0 / totalCost;
    }

    public void updateLocalBest(double fitness) {
        if (fitness > localBestFitness) {
            localBestFitness = fitness;
            bestVmAllocation = new ArrayList<>(vmAllocation);
            bestCloudletAllocation = new ArrayList<>(cloudletAllocation);
        }
    }

    public void updateGlobalBest(double fitness, List<Integer> vmAllocation, List<Integer> cloudletAllocation) {
        if (fitness > globalBestFitness) {
            globalBestFitness = fitness;
            bestVmAllocation = new ArrayList<>(vmAllocation);
            bestCloudletAllocation = new ArrayList<>(cloudletAllocation);
        }
    }

    public void updatePosition(List<Integer> newVmAllocation, List<Integer> newCloudletAllocation) {
        vmAllocation = new ArrayList<>(newVmAllocation);
        cloudletAllocation = new ArrayList<>(newCloudletAllocation);
    }

    public void updateVelocity(List<Integer> bestVmAllocation, List<Integer> bestCloudletAllocation, double inertiaWeight, double cognitiveWeight, double socialWeight) {
        for (int i = 0; i < vmAllocation.size(); i++) {
            int currentVm = vmAllocation.get(i);
            int currentCloudlet = cloudletAllocation.get(i);
            int bestVm = bestVmAllocation.get(i);
            int bestCloudlet = bestCloudletAllocation.get(i);

            double velocity = inertiaWeight * currentVm +
                    cognitiveWeight * (bestVm - currentVm) * Math.random() +
                    socialWeight * (bestCloudlet - currentCloudlet) * Math.random();

            vmAllocation.set(i, (int) Math.round(currentVm + velocity));
        }

        updatePosition(vmAllocation, cloudletAllocation);
    }

    public List<Integer> getBestVmAllocation() {
        return bestVmAllocation;
    }

    public List<Integer> getBestCloudletAllocation() {
        return bestCloudletAllocation;
    }
}
