package org.cloudbus.cloudsim.examples;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.cloudbus.cloudsim.Cloudlet;
import org.cloudbus.cloudsim.Vm;

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
        // The parameters below are based on the reference paper
        this.femaleImmigrateRate = 0;
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
        // Individual[] individuals = population.getIndividuals();
        Individual[] individuals = population.getIndividuals();

        for (int i = 0; i < population.size(); i++) {
            Individual individual = individuals[i];
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
                if (i < maleNomadSize) {
                    individual.setIsMale(true);
                }

                // If the individual is a pride lion and prideIndex % lionsInPrid <
                // maleLionsInPride, the individual is a male lion
                // Otherwise, the individual is a female lion
            } else {

                // Set the prideId starting from 1 for pride lions
                individual.setPrideId((i % prideNumber) + 1);
                int prideIndex = (int) (i / prideNumber);
                int lionsInPride = (int) (prideSize / prideNumber);
                int maleLionsInPride = (int) (lionsInPride * malePridePercentage);
                boolean isMale = prideIndex % lionsInPride < maleLionsInPride;
                individual.setIsMale(isMale);
            }
        }
        return population;
    }

    // 2) Apply Lion Algorithm
    public Population applyLionBehavior(Population population, int dataCenterIterator, int iteration) {
        // Apply behavior based on whether the individual is a nomad or part of a pride
        // Step 2.1) Hunt
        prideHunt(population, dataCenterIterator, iteration);

        // Step 2.2) Nomad Roam
        nomadRoam(population, dataCenterIterator, iteration);

        // Step 2.3) Nomad Attack
        nomadAttack(population);

        // Step 2.4) Pride Immigration (THERE'S SOMETHING WRONG IN HERE!!!!)
        prideImmigrate(population);

        // Step 2.5) Nomad Immigration
        nomadImmigrate(population);
        return population;
    }

    // 2.1) Apply Hunting on Pride
    public void prideHunt(Population population, int dataCenterIterator, int cloudletIteration) {
        // 2.1.2) Select Female Hunters on each prides
        Map<Integer, List<Individual>> huntersByPride = selectFemaleHunters(population);

        for (Individual individual : population.getIndividuals()) {
            // 2.1.3) Check if the individual belongs to a pride
            if (individual.getPrideId() > 0) {
                // 2.1.4) Select the female hunters
                List<Individual> hunters = huntersByPride.get(individual.getPrideId());

                // 2.1.5) Partition the hunters into three groups
                List<Individual> leftWing = new ArrayList<>();
                List<Individual> center = new ArrayList<>();
                List<Individual> rightWing = new ArrayList<>();
                partitionHunters(hunters, leftWing, center, rightWing, dataCenterIterator, cloudletIteration);

                // 2.1.6) Calculate the prey's position
                int[] prey = calculatePreyPosition(hunters, dataCenterIterator);

                // 2.1.7) Move each hunter towards the prey
                if (hunters == null) {
                    hunters = new ArrayList<>();
                }
                else {
                    for (Individual hunter : hunters) {
                        int[] oldPositions = Arrays.copyOf(hunter.getVmPositions(), hunter.getVmPositions().length);
                        int[] newPositions;
    
                        if (center.contains(hunter)) {
                            newPositions = calculateNewCenterHunterPosition(hunter, prey, dataCenterIterator);
                        } else {
                            newPositions = calculateNewWingHunterPosition(hunter, prey, dataCenterIterator);
                        }
    
                        double oldFitness = hunter.getFitness();
                        hunter.setVmPositions(newPositions);
                        double newFitness = calcFitness(hunter, dataCenterIterator, cloudletIteration);
                        // 2.1.8) The new position is not better, update the hunter's fitness to
                        // oldFitness
                        if (newFitness < oldFitness) {
                            hunter.setVmPositions(oldPositions);
                            hunter.setFitness(oldFitness);
                        }
                    }
                }  
            }
        }
    }

    // 2.1.2) Select Female Hunters on each prides
    public Map<Integer, List<Individual>> selectFemaleHunters(Population population) {
        Map<Integer, List<Individual>> individualsByPride = Arrays.stream(population.getIndividuals())
                .filter(lion -> !lion.getIsMale() && lion.getPrideId() > 0)
                .collect(Collectors.groupingBy(Individual::getPrideId));

        // 2.1.2.1) Initialize the map
        Map<Integer, List<Individual>> huntersByPride = new HashMap<>();

        // 2.1.2.2) For each pride, select hunters
        for (Map.Entry<Integer, List<Individual>> entry : individualsByPride.entrySet()) {
            List<Individual> females = new ArrayList<>(entry.getValue());
            Random rand = new Random();
            // 2.1.2.3) Randomly decide the number of hunters, at least one
            int numHunters = rand.nextInt(females.size()) + 1;

            // 2.1.2.4) Save the selected hunters in a list
            List<Individual> hunters = new ArrayList<>();
            for (int i = 0; i < numHunters; i++) {
                int index = rand.nextInt(females.size());
                hunters.add(females.get(index));
                // 2.1.2.5) Remove the selected lioness to avoid selecting it again
                females.remove(index);
            }

            huntersByPride.put(entry.getKey(), hunters);
        }

        return huntersByPride;
    }

    // 2.1.5) Partition the hunters into three groups
    public void partitionHunters(List<Individual> hunters, List<Individual> leftWing, List<Individual> center,
            List<Individual> rightWing, int dataCenterIterator, int cloudletIteration) {
        // 2.1.5.1) Sort the hunters by fitness in ascending order of x (or descending
        // order of fitness)
        if (hunters == null) {
            hunters = new ArrayList<>();
        }
        hunters.sort(Comparator.comparingDouble(Individual::getFitness));
        hunters.sort(Comparator.comparingDouble((Individual a) -> a.getFitness()));

        // 2.1.5.2) Partition the sorted hunters into three groups
        int size = hunters.size();
        // 2.1.5.3) One third of the hunters go to the center
        int centerSize = size / 3;
        // 2.1.5.4) The rest are evenly divided between the two wings
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

    // 2.1.6) Calculate the prey position
    public int[] calculatePreyPosition(List<Individual> hunters, int dataCenterIterator) {
        if (hunters == null || hunters.isEmpty() || hunters.get(0).getVmPositions().length == 0) {
            return new int[0];
        }
        else {
            int[] preyPositions = new int[hunters.get(0).getVmPositions().length];
            int[] totalPositions = new int[hunters.get(0).getVmPositions().length];
            int count = 0;

            int minGene = (dataCenterIterator) * 9;
            int maxGene = ((dataCenterIterator + 1) * 9) - 1;

            for (Individual hunter : hunters) {
                if (hunter.getFitness() > 0) {
                    int[] hunterPositions = hunter.getVmPositions();
                    for (int i = 0; i < hunterPositions.length; i++) {
                        totalPositions[i] += Math.abs(hunterPositions[i]);
                    }
                    count++;
                }
            }

            if (count != 0) {
                for (int i = 0; i < preyPositions.length; i++) {
                    preyPositions[i] = Math.abs(totalPositions[i] / count);
                    // Ensure preyPositions[i] is within the range of minGene and maxGene
                    // preyPositions[i] = Math.max(0, Math.min(62, preyPositions[i]));
                    preyPositions[i] = Math.max(minGene, Math.min(maxGene, preyPositions[i]));
                }
            }
            return preyPositions;
        }
    }

    // 2.1.7.1) Move the center hunter towards the prey
    public int[] calculateNewCenterHunterPosition(Individual hunter, int[] prey, int dataCenterIterator) {
        int[] hunterPositions = hunter.getVmPositions();
        int[] newPositions = new int[hunterPositions.length];
        Random rand = new Random();
        int minGene = (dataCenterIterator) * 9;
        int maxGene = (dataCenterIterator + 1) * 9 - 1;

        for (int i = 0; i < hunterPositions.length; i++) {
            if (hunterPositions[i] < prey[i]) {
                // int bound = Math.abs((prey[i] - hunterPositions[i]) + hunterPositions[i]);
                // if (bound <= 0) {

                //     newPositions[i] = rand.nextInt(1);
                // } else {
                //     newPositions[i] = rand.nextInt(bound);
                // }
                newPositions[i] = rand.nextInt(prey[i] - hunterPositions[i]) + hunterPositions[i];
            } else if (hunterPositions[i] > prey[i]) {
                // int bound = Math.abs((hunterPositions[i] - prey[i]) + prey[i]);
                // if (bound <= 0) {
                //     newPositions[i] = rand.nextInt(1);
                // } else {
                //     newPositions[i] = rand.nextInt(bound);
                // }
                newPositions[i] = rand.nextInt(hunterPositions[i] - prey[i]) + prey[i];
            }
            else {
                newPositions[i] = hunterPositions[i];
            }
            // newPositions[i] = Math.max(minGene, Math.min(maxGene, newPositions[i]));
        }

        return newPositions;
    }

    // 2.1.7.2) Move the wings hunter towards the prey
    public int[] calculateNewWingHunterPosition(Individual hunter, int[] prey, int dataCenterIterator) {
        int[] hunterPositions = hunter.getVmPositions();
        int[] newPositions = new int[hunterPositions.length];
        Random rand = new Random();
        int minGene = (dataCenterIterator) * 9;
        int maxGene = ((dataCenterIterator + 1) * 9) - 1;

        for (int i = 0; i < hunterPositions.length; i++) {
            int temp = Math.abs(2 * prey[i] - hunterPositions[i]);
            if (temp < prey[i]) {
                newPositions[i] = rand.nextInt((prey[i] - temp)) + temp;
            } 
//            else if (temp < prey[i]){
//                newPositions[i] = rand.nextInt((temp - prey[i])) + prey[i];
//            }
            else {
                newPositions[i] = hunterPositions[i];
            }
            // newPositions[i] = Math.max(minGene, Math.min(maxGene, newPositions[i]));
        }

        return newPositions;
    }

    // 2.2) Nomad Roam
    public void nomadRoam(Population population, int dataCenterIterator, int cloudletIteration) {
        Random rand = new Random();
        double pr = 0.5;
    
        for (Individual lion : population.getIndividuals()) {
            // 2.2.1) Check if the lion is a nomad
            if (lion.getPrideId() == 0) {
                // 2.2.2) Check if the lion should move
                if (rand.nextDouble() < pr) {
                    int[] oldPositions = lion.getVmPositions();
                    int[] newPositions = new int[oldPositions.length];
    
                    int minGene = (dataCenterIterator - 1) * 9;
                    int maxGene = ((dataCenterIterator + 1) * 9) - 1;
    
                    for (int i = 0; i < oldPositions.length; i++) {
                        // 2.2.3) Replace with your own logic for generating a random position
                        newPositions[i] = rand.nextInt(maxGene - minGene + 1) + minGene;
    
                        // Ensure the new position is within the range of minGene and maxGene
                        newPositions[i] = Math.max(minGene, Math.min(maxGene, newPositions[i]));
                    }
    
                    double oldFitness = lion.getFitness();
    
                    // 2.2.4) Calculate the fitness at the new position
                    lion.setVmPositions(newPositions);
                    double newFitness = calcFitness(lion, dataCenterIterator, cloudletIteration);
    
                    // 2.2.5) If the new position is not better, revert the lion's position and
                    // fitness
                    if (newFitness < oldFitness) {
                        lion.setVmPositions(oldPositions);
                        lion.setFitness(oldFitness);
                    }
                }
            }
        }
    }

    // 2.3) Nomads Attack the Pride Lions
    public void nomadAttack(Population population) {
        // 2.3.1) Method to get all nomad lions
        List<Individual> nomads = population.getNomads();
        Random rand = new Random();

        // 2.3.2) Random number of nomad attacks
        int numAttacks = rand.nextInt(nomads.size()) + 1;

        for (int i = 0; i < numAttacks; i++) {
            Individual nomadLion = nomads.get(rand.nextInt(nomads.size()));

            // 2.3.3) Ensure the nomad lion is male
            if (!nomadLion.getIsMale()) {
                continue;
            }

            // 2.3.4) Assuming a constant number of prides = 4
            int targetPrideId = rand.nextInt(prideNumber) + 1;
            // 2.3.5) Method to get the male lion with the lowest fitness from a pride
            Individual lowestFitnessPrideMaleLion = population.getLowestFitnessPrideMaleLion(targetPrideId);

            // 2.3.6) If the nomad lion has better fitness, it replaces the male pride lion
            if (nomadLion.getFitness() > lowestFitnessPrideMaleLion.getFitness()) {
                // 2.3.7) The nomad lion joins the pride
                nomadLion.setPrideId(targetPrideId);

                // 2.3.8) The replaced male pride lion becomes a nomad
                lowestFitnessPrideMaleLion.setPrideId(0);
                // 2.3.9) Ensure the lion remains male
                lowestFitnessPrideMaleLion.setIsMale(true);
            }
        }
    }

    // 2.4) Immigrate the female pride
    public void prideImmigrate(Population population) {
        for (int prideId = 1; prideId <= this.prideNumber; prideId++) {
            List<Individual> prideLions = population.getPrideLions(prideId);
            int numFemalesToMigrate = (int)(this.femaleImmigrateRate * prideLions.size());
            Random rand = new Random();

            for (int i = 0; i < numFemalesToMigrate; i++) {
                List<Individual> femaleLions = prideLions.stream().filter(individual -> !individual.getIsMale()).collect(Collectors.toList());
                if (!femaleLions.isEmpty()) {
                    int index = rand.nextInt(femaleLions.size());
                    Individual femaleLion = femaleLions.get(index);
                    femaleLion.setPrideId(0);  // Set nomadic status
                }
            }
        }
    }

// 2.5) Immigrate the female nomad
public void nomadImmigrate(Population population) {
    List<Individual> femaleNomads = Arrays.stream(population.getIndividuals())
        .filter(individual -> !individual.getIsMale() && individual.getPrideId() == 0)
        .sorted(Comparator.comparingDouble(Individual::getFitness))
        .collect(Collectors.toList());

    for (int prideId = 1; prideId <= prideNumber; prideId++) {
        long femaleCount = population.getPrideLions(prideId).stream().filter(individual -> !individual.getIsMale()).count();

        if (femaleCount < (int) ((1 - this.malePridePercentage) * population.getPrideLions(prideId).size())) {
            if (!femaleNomads.isEmpty()) {
                Individual highestFitnessFemaleNomad = femaleNomads.remove(0);
                highestFitnessFemaleNomad.setPrideId(prideId);
            }
        }
    }
}


    // // 2.4) Immigrate the female pride
    // public void prideImmigrate(Population population) {
    //     // 2.4.1) Iterate through each pride
    //     for (int prideId = 1; prideId <= this.prideNumber; prideId++) {
    //         List<Individual> prideLions = population.getPrideLions(prideId);
    //         // List<Individual> femaleLions = population.getPrideFemaleLions(prideId);
    //         Random rand = new Random();

    //         // 2.4.2) Separate female lions
    //         // femaleLions = femaleLions.stream()
    //         List<Individual> femaleLions = prideLions.stream().filter(individual -> !individual.getIsMale()).collect(Collectors.toList());

    //         // 2.4.3) Calculate the number of lions to migrate
    //         int numFemalesToMigrate = (int) (this.femaleImmigrateRate * femaleLions.size());
    //         // int numFemalesToMigrate = (int) (0 * femaleLions.size());
    //         // int numFemalesToMigrate = 4;

    //         // 2.4.4) Migrate the selected number of female lions
    //         for (int i = 0; i < numFemalesToMigrate; i++) {
    //             if (!femaleLions.isEmpty()) {
    //                 int index = rand.nextInt(femaleLions.size());
    //                 Individual femaleLion = femaleLions.get(index);
    //                 // 2.4.5) Setting prideId to 0 to indicate nomadic status
    //                 femaleLion.setPrideId(0);
    //                 femaleLions.remove(index);
    //             }
    //         }
    //     }
    // }

    // // 2.5) Immigrate the female nomad
    // public void nomadImmigrate(Population population) {
    //     // 2.5.1) Sort all female nomads by fitness value
    //     List<Individual> femaleNomads = Arrays.stream(population.getIndividuals())
    //         .filter(individual -> !individual.getIsMale() && individual.getPrideId() == 0)
    //         .sorted((individual1, individual2) -> Double.compare(individual2.getFitness(), individual1.getFitness()))
    //         .collect(Collectors.toList());

    //     // 2.5.2) Look for prides with female emptiness
    //     for (int prideId = 1; prideId <= prideNumber; prideId++) {
    //         List<Individual> prideLions = population.getPrideLions(prideId);
    //         int femaleCount = (int) prideLions.stream().filter(individual -> !individual.getIsMale()).count();
    //         int totalPrideSize = prideLions.size();

    //         // 2.5.3) Check if the pride follows the rules
    //         if ((totalPrideSize < (((1 - this.roamingPercentage) * this.populationSize)/this.prideNumber)) && (femaleCount < (int) ((1 - this.malePridePercentage) * totalPrideSize))) {
    //             // 2.5.4) Make the highest fitness female nomad enter the pride
    //             if (!femaleNomads.isEmpty()) {
    //                 Individual highestFitnessFemaleNomad = femaleNomads.get(femaleNomads.size() - 1);
    //                 highestFitnessFemaleNomad.setPrideId(prideId);
    //                 femaleNomads.remove(highestFitnessFemaleNomad);
    //             }
    //         }
    //     }
    // }

    public void evalPopulation(Population population, int dataCenterIterator, int cloudletIteration) {
        for (Individual individual : population.getIndividuals()) {
            double individualFitness = calcFitness(individual, dataCenterIterator, cloudletIteration);
            individual.setFitness(Math.abs(individualFitness));
        }
    }

    public double calcFitness(Individual individual, int dataCenterIterator, int cloudletIteration) {
        double totalExecutionTime = 0;
        double totalCost = 0;
        int iterator = 0;
        dataCenterIterator = dataCenterIterator - 1;

        for (int i = 0 + dataCenterIterator * 9 + cloudletIteration * 54; i < 9 + dataCenterIterator * 9
                + cloudletIteration * 54; i++) {
            int gene = Math.abs(individual.getGene(iterator));
            // int gene = individual.getGene(iterator);
            double mips = calculateMips(gene % 9);

            totalExecutionTime += cloudletList.get(i).getCloudletLength() / mips;
            totalCost += calculateCost(vmList.get(gene % 9), cloudletList.get(i));
            iterator++;
        }

        // Modify the fitness calculation to prioritize makespan optimization
        double makespanFitness = calculateMakespanFitness(totalExecutionTime);
        double costFitness = calculateCostFitness(totalCost);
        double fitness = makespanFitness + costFitness;

        individual.setFitness(Math.abs(fitness));
        return fitness;
    }

    public double calculateTotalFitness(List<Individual> individuals) {
        double totalFitness = 0;
        for (Individual individual : individuals) {
            totalFitness += individual.getFitness();
        }
        return totalFitness;
    }

    // public int applyLionBehaviorToGene(int gene, int dataCenterIterator) {
    // int minGene = dataCenterIterator * 9;
    // int maxGene = (dataCenterIterator + 1) * 9 - 1;

    // if (gene < minGene || gene > maxGene) {
    // return generateRandomNonEqualVmgene(gene, dataCenterIterator);
    // }

    // int squirrelJump = 1; // Adjust the squirrel jump value based on the problem
    // domain

    // // Apply squirrel movement by adding/subtracting squirrelJump to the current
    // // gene
    // int newGene = gene + squirrelJump;
    // if (newGene > maxGene) {
    // newGene = maxGene;
    // } else if (newGene < minGene) {
    // newGene = minGene;
    // }

    // return newGene;
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

    // public void updateGlobalBest(double fitness, List<Integer> vmPositions,
    // List<Integer> cloudletAllocation) {
    // if (fitness > globalBestFitness) {
    // globalBestFitness = fitness;
    // vmBestPositions = new ArrayList<>(vmPositions);
    // bestCloudletAllocation = new ArrayList<>(cloudletAllocation);
    // }
    // }

    // public void updatePosition(List<Integer> newvmPositions, List<Integer>
    // newCloudletAllocation) {
    // vmPositions = new ArrayList<>(newvmPositions);
    // cloudletAllocation = new ArrayList<>(newCloudletAllocation);
    // }

    // public void updateVelocity(List<Integer> vmBestPositions, List<Integer>
    // bestCloudletAllocation,
    // double inertiaWeight, double cognitiveWeight, double socialWeight) {
    // for (int i = 0; i < vmPositions.size(); i++) {
    // int currentVm = vmPositions.get(i);
    // int currentCloudlet = cloudletAllocation.get(i);
    // int bestVm = vmBestPositions.get(i);
    // int bestCloudlet = bestCloudletAllocation.get(i);

    // double velocity = inertiaWeight * currentVm +
    // cognitiveWeight * (bestVm - currentVm) * Math.random() +
    // socialWeight * (bestCloudlet - currentCloudlet) * Math.random();

    // vmPositions.set(i, (int) Math.round(currentVm + velocity));
    // }

    // updatePosition(vmPositions, cloudletAllocation);
    // }

    public List<Integer> getvmBestPositions() {
        return vmBestPositions;
    }

    public List<Integer> getBestCloudletAllocation() {
        return bestCloudletAllocation;
    }
}
