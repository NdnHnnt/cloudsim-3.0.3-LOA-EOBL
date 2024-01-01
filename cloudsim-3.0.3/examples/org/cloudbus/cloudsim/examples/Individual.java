package org.cloudbus.cloudsim.examples;

import java.util.Arrays;
import java.util.Random;

// chromosome = vmPositions
// chromosomeLength = cloudletLength
public class Individual {
    private int[] vmPositions;
    private double fitness = -1;
    private boolean isMale;
    private int prideId;

    public Individual(int[] vmPositions) {
        this.vmPositions = vmPositions;
    }

    public Individual(int cloudletLength, int dataCenterIterator) {
        this.vmPositions = new int[cloudletLength];
        dataCenterIterator = dataCenterIterator - 1;
        int max = 8 + 9 * dataCenterIterator;
        int min = 0 + 9 * dataCenterIterator;
        int range = max - min + 1;

        Random random = new Random();
        for (int gene = 0; gene < cloudletLength; gene++) {
            int rand = Math.abs(random.nextInt(range) + min);
            setGene(gene, rand);
        }
    }

    public int[] getVmPositions() {
        return this.vmPositions;
    }

    public int getCloudletLength() {
        return this.vmPositions.length;
    }

    public void setGene(int offset, int gene) {
        this.vmPositions[offset] = gene;
    }

    public int getGene(int offset) {
        return this.vmPositions[offset];
    }

    public void setFitness(double fitness) {
        this.fitness = fitness;
    }

    public double getFitness() {
        return this.fitness;
    }

    public void setIsMale(boolean isMale) {
        this.isMale = isMale;
    }

    public boolean getIsMale() {
        return this.isMale;
    }

    public void setPrideId(int prideId) {
        this.prideId = prideId;
    }

    public int getPrideId() {
        return this.prideId;
    }

    public String toString() {
        StringBuilder output = new StringBuilder();
        for (int gene = 0; gene < this.vmPositions.length; gene++) {
            output.append(this.vmPositions[gene]);
        }
        return output.toString();
    }

    public void setPosition(int index, int newPosition) {
        if (index >= 0 && index < vmPositions.length) {
            vmPositions[index] = newPosition;
        }
    }

    public void setVmPositions(int[] newPositions) {
        this.vmPositions = Arrays.copyOf(newPositions, newPositions.length);
    }
}

