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
    private double nomadPercentage;
    private double maleNomadPercentage;
    private double malePridePercentage;
    private double roamingPercentage;
    private double matingPercentage;
    private double femaleImmigrateRate;
    private Map<Integer, List<Individual>> huntersByPride;

    public LionOptimizationAlgorithm(int populationSize, List<Cloudlet> cloudletList,List<Vm> vmList) {
        this.huntersByPride = new HashMap<>();
        this.populationSize = populationSize;
        this.cloudletList = cloudletList;
        this.vmList = vmList;
        // The parameters below are based on the reference paper
        this.prideNumber = 4;
        this.roamingPercentage = 0.2;
        this.nomadPercentage = 0.2;
        this.malePridePercentage = 0.2;
        this.maleNomadPercentage = 1 - this.malePridePercentage;
        this.femaleImmigrateRate = 0.4;
        this.matingPercentage = 0.3;
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
    public Population applyLionBehavior(Population population, int dataCenterIterator, int cloudletIterator) {
        // Apply behavior based on whether the individual is a nomad or part of a pride
        // Step 2.1) HuntPride
        prideHunt(population, dataCenterIterator, cloudletIterator);

        // Step 2.2) Pride Roam
        prideRoam(population, dataCenterIterator, cloudletIterator);

        // Step 2.3) Mating on Pride
        prideMating(population, dataCenterIterator, cloudletIterator);

        // Step 2.4) Nomad Roam
        nomadRoam(population, dataCenterIterator, cloudletIterator);

        // Step 2.5) Mating on Nomad
        nomadMating(population, dataCenterIterator, cloudletIterator);

        // Step 2.6) Nomad Attack
        nomadAttack(population);

        // Step 2.7) Pride Immigration
         prideImmigrate(population);

        // Step 2.8) Nomad Immigration
         nomadImmigrate(population);

        // Step 2.9) Apply Equilibrium
         lionEquilibrium(population);
        
        return population;
    }

    // 2.1) Apply Hunting on Pride
    public void prideHunt(Population population, int dataCenterIterator, int cloudletIteration) {
        // 2.1.2) Select Female Hunters on each pride
        Map<Integer, List<Individual>> huntersByPride = selectFemaleHunters(population);
    
        // Iterate over each pride
        for (int prideId = 1; prideId <= this.prideNumber; prideId++) {
            // 2.1.4) Select the female hunters for the current pride
            List<Individual> hunters = huntersByPride.get(prideId);
    
            // Skip the process if there are no hunters in the current pride
            if (hunters == null || hunters.isEmpty()) {
                continue;
            }
    
            // 2.1.5) Partition the hunters into three groups
            List<Individual> leftWing = new ArrayList<>();
            List<Individual> center = new ArrayList<>();
            List<Individual> rightWing = new ArrayList<>();
            partitionHunters(hunters, leftWing, center, rightWing, dataCenterIterator, cloudletIteration);
    
            // 2.1.6) Calculate the prey's position
            int[] prey = calculatePreyPosition(hunters, dataCenterIterator);
    
            // 2.1.7) Move each hunter towards the prey
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
    
                // 2.1.8) If the new position is not better, revert to old position and fitness
                if (newFitness < oldFitness) {
                    hunter.setVmPositions(oldPositions);
                    hunter.setFitness(oldFitness);
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
            int numHunters = rand.nextInt(females.size());

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
        hunters.sort(Comparator.comparingDouble((Individual a) -> a.getFitness()).reversed());

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

            // int minGene = (dataCenterIterator) * 9;
            // int maxGene = ((dataCenterIterator + 1) * 9) - 1;

            for (Individual hunter : hunters) {
                if (hunter.getFitness() > 0) {
                    int[] hunterPositions = hunter.getVmPositions();
                    for (int i = 0; i < hunterPositions.length; i++) {
                        totalPositions[i] += (hunterPositions[i]);
                    }
                    count++;
                }
            }

            if (count != 0) {
                for (int i = 0; i < preyPositions.length; i++) {
                    preyPositions[i] = (totalPositions[i] / count);
                    // Ensure preyPositions[i] is within the range of minGene and maxGene
                    // preyPositions[i] = Math.max(0, Math.min(62, preyPositions[i]));
                    // preyPositions[i] = Math.max(minGene, Math.min(maxGene, preyPositions[i]));
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

        for (int i = 0; i < hunterPositions.length; i++) {
            if (hunterPositions[i] < prey[i]) {
                newPositions[i] = rand.nextInt(prey[i] - hunterPositions[i]) + hunterPositions[i];
            } else if (hunterPositions[i] > prey[i]) {
                newPositions[i] = rand.nextInt(hunterPositions[i] - prey[i]) + prey[i];
            }
            else {
                newPositions[i] = hunterPositions[i];
            }
        }
        return newPositions;
    }

    // 2.1.7.2) Move the wings hunter towards the prey
    public int[] calculateNewWingHunterPosition(Individual hunter, int[] prey, int dataCenterIterator) {
        int[] hunterPositions = hunter.getVmPositions();
        int[] newPositions = new int[hunterPositions.length];
        Random rand = new Random();

        for (int i = 0; i < hunterPositions.length; i++) {
            int temp = Math.abs((2 * prey[i]) - hunterPositions[i]);
            if (temp < prey[i]) {
                newPositions[i] = rand.nextInt((prey[i] - temp)) + temp;
            } 
            // else if (temp > prey[i]){
            //     newPositions[i] = rand.nextInt((temp - prey[i])) + prey[i];
            // }
            else {
                newPositions[i] = hunterPositions[i];
            }
            // newPositions[i] = Math.max(minGene, Math.min(maxGene, newPositions[i]));
        }

        return newPositions;
    }

    // 2.2) Pride Roam
    public void prideRoam(Population population, int dataCenterIterator, int cloudletIteration) {
        Random rand = new Random();
    
        for (Individual lion : population.getIndividuals()) {
            // 2.2.1) Check if the lion is a nomad
            if (lion.getPrideId() > 0 && lion.getIsMale()) {
                // 2.2.2) Check if the lion should move
                if (rand.nextDouble() <= this.roamingPercentage) {
                    int[] oldPositions = lion.getVmPositions();
                    int[] newPositions = new int[oldPositions.length];
    
                    int minGene = (dataCenterIterator - 1) * 9;
                    int maxGene = ((dataCenterIterator) * 9) - 1;
    
                    for (int i = 0; i < oldPositions.length; i++) {
                        // 2.2.3) Replace with your own logic for generating a random position
                        newPositions[i] = rand.nextInt(maxGene - minGene) + minGene;
                    }
    
                    double oldFitness = lion.getFitness();
    
                    // 2.4.4) Calculate the fitness at the new position
                    lion.setVmPositions(newPositions);
                    double newFitness = calcFitness(lion, dataCenterIterator, cloudletIteration);
    
                    // 2.4.5) If the new position is not better, revert the lion's position and fitness
                    if (newFitness < oldFitness) {
                        lion.setVmPositions(oldPositions);
                        lion.setFitness(oldFitness);
                    }
                }
            }
        }
    }

    // 2.3) Pride Mating
    public void prideMating(Population population, int dataCenterIterator, int cloudletIteration) {

        // 2.3.1) Get the prides from the population
        List<Individual[]> prides = Arrays.stream(population.getIndividuals())
            .filter(individual -> individual.getPrideId() > 0)
            .collect(Collectors.groupingBy(Individual::getPrideId))
            .values()
            .stream()
            .map(list -> list.toArray(new Individual[0]))
            .collect(Collectors.toList());   
        
        // 2.3.2) Apply mating operation to the prides
            for (Individual[] pride : prides) {
                // 2.3.3) Ensure there are at least two individuals in the pride
                if (pride.length < 2) {
                    continue;
                }

                // 2.3.4) Separate males and females in the pride
                Individual[] males = Arrays.stream(pride)
                .filter(individual -> individual.getIsMale())
                .toArray(Individual[]::new);
                Individual[] females = Arrays.stream(pride)
                .filter(individual -> !individual.getIsMale())
                .toArray(Individual[]::new);

                // 2.3.5) Ensure there is at least one male and one female in the pride
                if (males.length == 0 || females.length == 0) {
                    continue;
                }

                Random rand = new Random();
                if (rand.nextDouble() <= this.matingPercentage) { 
                    // 2.3.6) Select one male and one female for mating
                    Individual parent1 = males[new Random().nextInt(males.length)];
                    Individual parent2 = females[new Random().nextInt(females.length)];

                    // 2.3.7) Create two new individuals
                    Individual offspring1 = new Individual(parent1.getCloudletLength(), dataCenterIterator);
                    Individual offspring2 = new Individual(parent1.getCloudletLength(), dataCenterIterator);

                    // 2.3.8) Determine a crossover point
                    int crossoverPoint = new Random().nextInt(parent1.getCloudletLength());

                    // 2.3.9) Perform crossover
                    for (int j = 0; j < parent1.getCloudletLength(); j++) {
                        if (j < crossoverPoint) {
                            offspring1.setGene(j, parent1.getGene(j));
                            offspring2.setGene(j, parent2.getGene(j));
                        } else {
                            offspring1.setGene(j, parent1.getGene(j));
                            offspring2.setGene(j, parent2.getGene(j));
                        }
                    }

                    // Set the gender of the offspring and their pride
                    offspring1.setIsMale(true);
                    offspring2.setIsMale(false);
                    offspring1.setPrideId(parent1.getPrideId());
                    offspring2.setPrideId(parent1.getPrideId());
                    }
                }     
        }

       
    // 2.4) Nomad Roam
    public void nomadRoam(Population population, int dataCenterIterator, int cloudletIteration) {
        Random rand = new Random();
    
        for (Individual lion : population.getIndividuals()) {
            // 2.4.1) Check if the lion is a nomad
            if (lion.getPrideId() == 0) {
                // 2.4.2) Check if the lion should move
                if (rand.nextDouble() <= this.roamingPercentage) {
                    int[] oldPositions = lion.getVmPositions();
                    int[] newPositions = new int[oldPositions.length];
    
                    int minGene = (dataCenterIterator - 1) * 9;
                    int maxGene = ((dataCenterIterator) * 9) - 1;
    
                    for (int i = 0; i < oldPositions.length; i++) {
                        // 2.4.3) Replace with your own logic for generating a random position
                        newPositions[i] = rand.nextInt(maxGene - minGene) + minGene;
                    }
    
                    double oldFitness = lion.getFitness();
    
                    // 2.4.4) Calculate the fitness at the new position
                    lion.setVmPositions(newPositions);
                    double newFitness = calcFitness(lion, dataCenterIterator, cloudletIteration);
    
                    // 2.4.5) If the new position is not better, revert the lion's position and
                    // fitness
                    if (newFitness < oldFitness) {
                        lion.setVmPositions(oldPositions);
                        lion.setFitness(oldFitness);
                    }
                }
            }
        }
    }

    // 2.5) Mating on Nomad
    public void nomadMating(Population population, int dataCenterIterator, int cloudletIteration) {
        // 2.5.1) Get the nomad lions from the population
        Individual[] nomadLions = Arrays.stream(population.getIndividuals())
        .filter(individual -> individual.getPrideId() == 0)
        .toArray(Individual[]::new);

        // 2.5.2) Ensure there are at least two nomad lions
        if (nomadLions.length < 2) {
            return;
        }

        // 2.5.3) Separate males and females among the nomad lions
        Individual[] males = Arrays.stream(nomadLions)
        .filter(Individual::getIsMale)
        .toArray(Individual[]::new);
        Individual[] females = Arrays.stream(nomadLions)
        .filter(individual -> !individual.getIsMale())
        .toArray(Individual[]::new);

        // 2.5.4) Ensure there is at least one male and one female among the nomad lions
        if (males.length == 0 || females.length == 0) {
            return;
        }

        Random rand = new Random();
        if (rand.nextDouble() <= this.matingPercentage) { 

            // 2.5.5) Select one male and one female for mating
            Individual parent1 = males[new Random().nextInt(males.length)];
            Individual parent2 = females[new Random().nextInt(females.length)];

            // 2.5.6) Create two new individuals
            Individual offspring1 = new Individual(parent1.getCloudletLength(), dataCenterIterator);
            Individual offspring2 = new Individual(parent1.getCloudletLength(), dataCenterIterator);

            // 2.5.7) Perform crossover
            int crossoverPoint = new Random().nextInt(parent1.getCloudletLength());
            for (int j = 0; j < parent1.getCloudletLength(); j++) {
                if (j < crossoverPoint) {
                    offspring1.setGene(j, parent1.getGene(j));
                    offspring2.setGene(j, parent2.getGene(j));
                } else {
                    offspring1.setGene(j, parent1.getGene(j));
                    offspring2.setGene(j, parent2.getGene(j));
                }
            }

            // 2.5.8) Perform crossover
            offspring1.setIsMale(true);
            offspring2.setIsMale(false);
            offspring1.setPrideId(parent1.getPrideId());
            offspring2.setPrideId(parent1.getPrideId());       
        } 
    }
    
    // 2.6) Nomads Attack the Pride Lions
    public void nomadAttack(Population population) {
        // 2.6.1) Method to get all nomad lions
        List<Individual> nomads = population.getNomads();
        Random rand = new Random();

        // 2.6.2) Random number of nomad attacks
        int numAttacks = rand.nextInt(nomads.size()) + 1;

        for (int i = 0; i < numAttacks; i++) {
            Individual nomadLion = nomads.get(rand.nextInt(nomads.size()));

            // 2.6.3) Ensure the nomad lion is male
            if (!nomadLion.getIsMale()) {
                continue;
            }

            // 2.6.4) Assuming a constant number of prides = 4
            int targetPrideId = rand.nextInt(prideNumber) + 1;
            // 2.6.5) Method to get the male lion with the lowest fitness from a pride
            Individual lowestFitnessPrideMaleLion = population.getLowestFitnessPrideMaleLion(targetPrideId);

            // 2.6.6) If the nomad lion has better fitness, it replaces the male pride lion
            if (nomadLion.getFitness() > lowestFitnessPrideMaleLion.getFitness()) {
                // 2.6.7) The nomad lion joins the pride
                nomadLion.setPrideId(targetPrideId);

                // 2.6.8) The replaced male pride lion becomes a nomad
                lowestFitnessPrideMaleLion.setPrideId(0);
                // 2.6.9) Ensure the lion remains male
                lowestFitnessPrideMaleLion.setIsMale(true);
            }
        }
    }

    // 2.7) Immigrate the female pride
    public void prideImmigrate(Population population) {
        for (int prideId = 1; prideId <= this.prideNumber; prideId++) {
            List<Individual> prideLions = population.getPrideLions(prideId);

            List<Individual> femaleLions = prideLions.stream()
                .filter(individual -> !individual.getIsMale())
                .collect(Collectors.toList());

            int numFemalesToMigrate = (int)(this.femaleImmigrateRate * femaleLions.size());
            int femalePrideCount = (int) population.getPrideLions(prideId).stream().filter(individual -> !individual.getIsMale()).count();
            int femalePrideIdeal = (int) (((this.populationSize * (1 - this.nomadPercentage))/this.prideNumber) * (1 - this.malePridePercentage));
            if (femalePrideCount > femalePrideIdeal) {
                numFemalesToMigrate += (femalePrideCount - femalePrideIdeal);
            }
            Random rand = new Random();

            List<Individual> malePride = Arrays.stream(population.getIndividuals())
                .filter(individual -> individual.getIsMale() && individual.getPrideId() > 0)
                .sorted(Comparator.comparingDouble(Individual::getFitness).reversed())
                .collect(Collectors.toList());

            for (int i = 0; i < numFemalesToMigrate; i++) {
                if (!femaleLions.isEmpty()) {
                    int index = rand.nextInt(femaleLions.size());
                    Individual femaleLion = femaleLions.get(index);
                    femaleLion.setPrideId(0);  // Set nomadic status
                }
            }

            int malePrideCount = (int) population.getPrideLions(prideId).stream().filter(individual -> individual.getIsMale()).count();
            if (malePrideCount > (int) (this.malePridePercentage * (this.populationSize * (1 - this.nomadPercentage)))){
                if (!malePride.isEmpty()) {
                    Individual lowestFitnessMalePride = malePride.get(0);
                    lowestFitnessMalePride.setPrideId(0);
                }
            }
        }
    }

    // 2.8) Immigrate the female nomad
    public void nomadImmigrate(Population population) {
        List<Individual> femaleNomads = Arrays.stream(population.getIndividuals())
            .filter(individual -> !individual.getIsMale() && individual.getPrideId() == 0)
            .sorted(Comparator.comparingDouble(Individual::getFitness))
            .collect(Collectors.toList());
        
        for (int prideId = 1; prideId <= prideNumber; prideId++) {
            int femalePrideCount = (int) population.getPrideLions(prideId).stream().filter(individual -> !individual.getIsMale()).count();
            
            if (femalePrideCount < (int) ((1 - this.malePridePercentage) * population.getPrideLions(prideId).size())) {
                if (!femaleNomads.isEmpty()) {
                	Individual highestFitnessFemaleNomad = femaleNomads.get(0);
                    highestFitnessFemaleNomad.setPrideId(prideId);
                }
            }
        }
    }

    // 2.9) Equilibrium Check
    public void lionEquilibrium(Population population) {
        List<Individual> femaleNomads = population.getNomads()
            .stream()
            .filter(individual -> !individual.getIsMale())
            .sorted(Comparator.comparingDouble(Individual::getFitness).reversed())
            .collect(Collectors.toList());
        
        List<Individual> maleNomads = population.getNomads()
            .stream()
            .filter(individual -> individual.getIsMale())
            .sorted(Comparator.comparingDouble(Individual::getFitness).reversed())
            .collect(Collectors.toList());
        
        int currentNomads = population.getNomads().size();
        int idealNomads = (int) (this.nomadPercentage * this.populationSize);       

        if (currentNomads > idealNomads){
            if (femaleNomads.size() > (this.malePridePercentage * idealNomads)){
                int numFemalesToRemove = (int) (femaleNomads.size() - (this.malePridePercentage * idealNomads));
                for (int i = 0; i < numFemalesToRemove; i++) {
                    Individual femaleNomad = femaleNomads.get(0);
                    femaleNomad.setPrideId(-1);
                    femaleNomad = femaleNomads.remove(0);
                }
            }
            if (maleNomads.size() > ((1 - this.malePridePercentage) * idealNomads)){
                int numMalesToRemove = (int)(maleNomads.size() - ((1 - this.malePridePercentage) * idealNomads));
                for (int i = 0; i < numMalesToRemove; i++) {
                    Individual maleNomad = maleNomads.get(0);
                    maleNomad.setPrideId(-1);
                    maleNomad = maleNomads.remove(0);
                }
            }
        }
    }

    // Step 2.10) Apply Elite Opposition-Based Learning
    public Individual eliteOppositionBasedLearning(Individual bestLion, int dataCenterIterator, int cloudletIterator){   
        int minGene = (dataCenterIterator - 1) * 9;
        int maxGene = ((dataCenterIterator + 1) * 9) - 1;
    
        // Save the initial fitness of the best lion
        double initialFitness = bestLion.getFitness();
    
        // Apply opposition-based learning
        for (int gene = 0; gene < bestLion.getCloudletLength(); gene++) {
            int opposite = maxGene - (bestLion.getGene(gene) - minGene);
            bestLion.setGene(gene, opposite);
        }
    
        // Recalculate the fitness of the best lion
        double newFitness = calcFitness(bestLion, dataCenterIterator, cloudletIterator);
    
        // Compare the fitness of the bSest lion before and after opposition-based learning
        if (newFitness < initialFitness) {
            // If the fitness has decreased, revert the changes
            for (int gene = 0; gene < bestLion.getCloudletLength(); gene++) {
                int original = minGene + maxGene - bestLion.getGene(gene);
                bestLion.setGene(gene, original);
            }
        } else {
            // If the fitness has increased, update the fitness of the best lion
            bestLion.setFitness(newFitness);
        }
    
        return bestLion;
    }

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
}
