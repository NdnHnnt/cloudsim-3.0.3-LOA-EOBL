package org.cloudbus.cloudsim.examples;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.cloudbus.cloudsim.Cloudlet;
import org.cloudbus.cloudsim.Vm;
import org.cloudbus.cloudsim.examples.Population;

public class LionOptimizationAlgorithm {
    private int prideNumber;
    private int populationSize;
    private List<Cloudlet> cloudletList;
    private List<Vm> vmList;
    private List<Integer> vmPositions;
    private List<Integer> vmBestPositions;
    private List<Integer> cloudletAllocation;
    private List<Integer> bestCloudletAllocation;
    private double localBestFitness;
    private double globalBestFitness;
    private double mutationRate;
    private double nomadPercentage;
    private double maleNomadPercentage;
    private double malePridePercentage;
    private double roamingPercentage;
    private double matingPercentage;
    private double femaleImmigrateRate;
    private Map<Integer, List<Individual>> huntersByPride;

    public LionOptimizationAlgorithm(int populationSize, double mutationRate, List<Cloudlet> cloudletList,
            List<Vm> vmList) {
        this.huntersByPride = new HashMap<>();
        this.populationSize = populationSize;
        this.mutationRate = mutationRate;
        this.cloudletList = cloudletList;
        this.vmList = vmList;
        this.localBestFitness = 0;
        this.globalBestFitness = 0;
        // The parameters below are based on the reference paper
        this.prideNumber = 4;
        this.roamingPercentage = 0.2;
        this.matingPercentage = 0.3;
        this.nomadPercentage = 0.2;
        this.malePridePercentage = 0.2;
        this.maleNomadPercentage = 1 - this.malePridePercentage;
    }

    // 1) Initialize the Population of Lions
    public Population initPopulation(int cloudletNumber, int dataCenterIterator) {
        Population population = new Population(this.populationSize, cloudletNumber, dataCenterIterator);

        // Calculate the number of nomad lions and the number of male and female lions
        // in both pride and nomad
        int nomadSize = (int) (this.nomadPercentage * this.populationSize);
        int prideSize = this.populationSize - nomadSize;
        int maleNomadSize = (int) (this.maleNomadPercentage * nomadSize);

        // Assign lion types (pride or nomad - male or female) to the individuals
        List<Individual> individuals = population.getIndividuals();

        for (int i = 0; i < this.populationSize; i++) {
            Individual individual = individuals.get(i);
            // If the individual's index is less than nomadSize, the individual is a nomad
            // lion
            // Otherwise, the individual is a pride lion
            boolean isNomad = i < nomadSize;

            // If the individual is a nomad lion and its index is less than maleNomadSize,
            // the individual is a male lion
            // Otherwise, the individual is a female lion
            if (isNomad) {

                // Set the prideId to 0 for nomad lions
                individual.setPrideId(0);
                boolean isMale = i < maleNomadSize;
                individual.setIsMale(isMale);

                // If the individual is a pride lion and prideIndex % lionsInPrid <
                // maleLionsInPride, the individual is a male lion
                // Otherwise, the individual is a female lion
            } else {

                // Set the prideId starting from 1 for pride lions
                individual.setPrideId((i % 4) + 1);
                int prideIndex = i / 4;
                int lionsInPride = (prideSize / 4);
                int maleLionsInPride = (int) (lionsInPride * malePridePercentage);
                boolean isMale = prideIndex % lionsInPride < maleLionsInPride;
                individual.setIsMale(isMale);
            }
        }
        return population;
    }

    public void prideHunt(Population population) {
        Map<Integer, List<Individual>> huntersByPride = selectFemaleHunters(population);

        for (Individual individual : population.getIndividuals()) {
            // Check if the individual belongs to a pride
            if (individual.getPrideId() > 0) { 
                // Select the female hunters
                List<Individual> hunters = huntersByPride.get(individual.getPrideId());

                // Partition the hunters into three groups
                List<Individual> leftWing = new ArrayList<>();
                List<Individual> center = new ArrayList<>();
                List<Individual> rightWing = new ArrayList<>();
                partitionHunters(hunters, leftWing, center, rightWing);

                // Calculate the prey's position
                double prey = calculatePreyPosition(hunters);

                // Move each hunter towards the prey
                for (Individual hunter : hunters) {
                    if (center.contains(hunter)) {
                        moveCenterHunter(hunter, prey);
                    } else {
                        moveWingHunter(hunter, prey);
                    }

                    // If the new position is better, update the prey's position
                    double newFitness = calcFitness(hunter, individual.getDataCenterIterator(), individual.getCloudletIteration());
                    if (newFitness > hunter.getFitness()) {
                        prey = updatePreyPosition(prey, hunter);
                        hunter.setFitness(newFitness);
                    }
                }
            }
        }
    }

    public Map<Integer, List<Individual>> selectFemaleHunters(Population population) {
        Map<Integer, List<Individual>> individualsByPride = Arrays.stream(population.getIndividuals())
                .filter(lion -> !lion.getIsMale() && lion.getPrideId() > 0)
                .collect(Collectors.groupingBy(Individual::getPrideId));

        // For each pride, select hunters
        for (Map.Entry<Integer, List<Individual>> entry : individualsByPride.entrySet()) {
            List<Individual> females = new ArrayList<>(entry.getValue());
            int numHunters = rand.nextInt(females.size()) + 1; // Randomly decide the number of hunters, at least one

            List<Individual> hunters = new ArrayList<>();
            for (int i = 0; i < numHunters; i++) {
                int index = rand.nextInt(females.size());
                hunters.add(females.get(index));
                females.remove(index); // Remove the selected lioness to avoid selecting it again
            }

            this.huntersByPride.put(entry.getKey(), hunters);
        }

        return huntersByPride;
    }

    public void partitionHunters(List<Individual> hunters, List<Individual> leftWing, List<Individual> center,
            List<Individual> rightWing) {
        // Sort the hunters by fitness in descending order
        hunters.sort((a, b) -> Double.compare(calcFitness(b), calcFitness(a)));

        // Partition the sorted hunters into three groups
        int size = hunters.size();
        // One third of the hunters go to the center
        int centerSize = size / 3;
        // The rest are evenly divided between the two wings
        int wingSize = (size - centerSize) / 2;

        for (int i = 0; i < size; i++) {
            if (i < centerSize) {
                center.add(hunters.get(i));
            } else if (i < centerSize + wingSize) {
                leftWing.add(hunters.get(i));
            } else {
                rightWing.add(hunters.get(i));
            }
        }
    }

    public double calcFitness(Individual individual, int dataCenterIterator, int cloudletIteration) {
        double totalExecutionTime = 0;
        double totalCost = 0;
        int iterator = 0;
        dataCenterIterator = dataCenterIterator - 1;

        for (int i = 0 + dataCenterIterator * 9 + cloudletIteration * 54; i < 9 + dataCenterIterator * 9
                + cloudletIteration * 54; i++) {
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

    public double calculateTotalFitness(List<Individual> individuals) {
        double totalFitness = 0;
        for (Individual individual : individuals) {
            totalFitness += individual.getFitness();
        }
        return totalFitness;
    }

    public double calculatePreyPosition(List<Individual> hunters) {
        double totalPosition = 0;
        for (Individual hunter : hunters) {
            totalPosition += hunter.getPosition();
        }
        return totalPosition / hunters.size();
    }

    public void moveHunter(Individual hunter, double prey, boolean isCenter) {
        if (isCenter) {
            if (hunter.getPosition() < prey) {
                hunter.setPosition(rand(hunter.getPosition(), prey));
            } else {
                hunter.setPosition(rand(prey, hunter.getPosition()));
            }
        } else {
            double value = 2 * prey - hunter.getPosition();
            if (value < prey) {
                hunter.setPosition(rand(value, prey));
            } else {
                hunter.setPosition(rand(prey, value));
            }
        }
    }

    public double updatePreyPosition(double prey, Individual hunter, double improvementPercentage) {
        return prey + Math.random() * 0.1 * improvementPercentage * Math.PI * (prey - hunter.getPosition());
    }

    public void evalPopulation(Population population, int dataCenterIterator, int cloudletIteration) {
        for (Individual individual : population.getIndividuals()) {
            double individualFitness = calcFitness(individual, dataCenterIterator, cloudletIteration);
            individual.setFitness(individualFitness);
        }
    }

    // public double calcFitness(Individual individual, int dataCenterIterator, int cloudletIteration) {
    //     double totalExecutionTime = 0;
    //     double totalCost = 0;
    //     int iterator = 0;
    //     dataCenterIterator = dataCenterIterator - 1;

    //     for (int i = 0 + dataCenterIterator * 9 + cloudletIteration * 54; i < 9 + dataCenterIterator * 9
    //             + cloudletIteration * 54; i++) {
    //         int gene = individual.getGene(iterator);
    //         double mips = calculateMips(gene % 9);

    //         totalExecutionTime += cloudletList.get(i).getCloudletLength() / mips;
    //         totalCost += calculateCost(vmList.get(gene % 9), cloudletList.get(i));
    //         iterator++;
    //     }

    //     // Modify the fitness calculation to prioritize makespan optimization
    //     double makespanFitness = calculateMakespanFitness(totalExecutionTime);
    //     double costFitness = calculateCostFitness(totalCost);
    //     double fitness = makespanFitness + costFitness;

    //     individual.setFitness(fitness);
    //     return fitness;
    // }

    // public void mutatePopulation(Population population, int dataCenterIterator) {
    //     Individual[] individuals = population.getIndividuals();

    //     for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
    //         Individual individual = individuals[populationIndex];

    //         for (int geneIndex = 0; geneIndex < individual.getChromosomeLength(); geneIndex++) {
    //             if (Math.random() < this.mutationRate) {
    //                 // Randomly select a new VM gene that is not equal to the current gene
    //                 int currentGene = individual.getGene(geneIndex);
    //                 int newGene = generateRandomNonEqualVmgene(currentGene, dataCenterIterator);
    //                 individual.setGene(geneIndex, newGene);
    //             }
    //         }
    //     }
    // }

    // public Population applyLionBehavior(Population population, int dataCenterIterator) {
    //     Individual[] individuals = population.getIndividuals();

    //     for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
    //         Individual individual = individuals[populationIndex];

    //         for (int geneIndex = 0; geneIndex < individual.getChromosomeLength(); geneIndex++) {
    //             if (Math.random() < this.mutationRate) {
    //                 int currentGene = individual.getGene(geneIndex);
    //                 int newGene = applyLionBehaviorToGene(currentGene, dataCenterIterator);
    //                 individual.setGene(geneIndex, newGene);
    //             }
    //         }
    //     }

    //     return population;
    // }

    // public int generateRandomNonEqualVmgene(int currentGene, int dataCenterIterator) {
    //     int minGene = dataCenterIterator * 9;
    //     int maxGene = (dataCenterIterator + 1) * 9 - 1;
    //     int newGene;

    //     // Generate a random VM gene until it is not equal to the current gene
    //     do {
    //         newGene = new Random().nextInt(maxGene - minGene + 1) + minGene;
    //     } while (newGene == currentGene);

    //     return newGene;
    // }

    // public int applyLionBehaviorToGene(int gene, int dataCenterIterator) {
    //     int minGene = dataCenterIterator * 9;
    //     int maxGene = (dataCenterIterator + 1) * 9 - 1;

    //     if (gene < minGene || gene > maxGene) {
    //         return generateRandomNonEqualVmgene(gene, dataCenterIterator);
    //     }

    //     int squirrelJump = 1; // Adjust the squirrel jump value based on the problem domain

    //     // Apply squirrel movement by adding/subtracting squirrelJump to the current
    //     // gene
    //     int newGene = gene + squirrelJump;
    //     if (newGene > maxGene) {
    //         newGene = maxGene;
    //     } else if (newGene < minGene) {
    //         newGene = minGene;
    //     }

    //     return newGene;
    // }

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
            vmBestPositions = new ArrayList<>(vmPositions);
            bestCloudletAllocation = new ArrayList<>(cloudletAllocation);
        }
    }

    public void updateGlobalBest(double fitness, List<Integer> vmPositions, List<Integer> cloudletAllocation) {
        if (fitness > globalBestFitness) {
            globalBestFitness = fitness;
            vmBestPositions = new ArrayList<>(vmPositions);
            bestCloudletAllocation = new ArrayList<>(cloudletAllocation);
        }
    }

    public void updatePosition(List<Integer> newvmPositions, List<Integer> newCloudletAllocation) {
        vmPositions = new ArrayList<>(newvmPositions);
        cloudletAllocation = new ArrayList<>(newCloudletAllocation);
    }

    public void updateVelocity(List<Integer> vmBestPositions, List<Integer> bestCloudletAllocation,
            double inertiaWeight, double cognitiveWeight, double socialWeight) {
        for (int i = 0; i < vmPositions.size(); i++) {
            int currentVm = vmPositions.get(i);
            int currentCloudlet = cloudletAllocation.get(i);
            int bestVm = vmBestPositions.get(i);
            int bestCloudlet = bestCloudletAllocation.get(i);

            double velocity = inertiaWeight * currentVm +
                    cognitiveWeight * (bestVm - currentVm) * Math.random() +
                    socialWeight * (bestCloudlet - currentCloudlet) * Math.random();

            vmPositions.set(i, (int) Math.round(currentVm + velocity));
        }

        updatePosition(vmPositions, cloudletAllocation);
    }

    public List<Integer> getvmBestPositions() {
        return vmBestPositions;
    }

    public List<Integer> getBestCloudletAllocation() {
        return bestCloudletAllocation;
    }
}
