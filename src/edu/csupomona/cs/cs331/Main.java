package edu.csupomona.cs.cs331;

/**
 * Objective:
 * Use and compare run times of three different methods
 * of matrix multiplication
 * 
 * Prof:	Dr. Young
 * Class:	CS 331
 * @author	Andrew Armaneous
 */

public class Main {

	public static void main(String[] args) {
		Matrix matrix;
		int size = 2;					// Initial size of matrices
		long startTime, endTime;		// For determining run times
		int round = 1;					// Current round of matrix size increase
		
		while(true){
			// Initialize matrices using size
			matrix = new Matrix(size);
			System.out.println();
			System.out.println("---------------------------------------");
			System.out.println("Round number: " + round);
			System.out.println("---------------------------------------");
			
			// Compute using 'naive', or classical, matrix multiplication
			startTime = System.nanoTime();
			matrix.naiveMultiply();
			endTime = System.nanoTime();
			// Print computation times in milliseconds
			System.out.println("Naive Multiplication: " + (endTime - startTime));
			
			// Compute using Divide and Conquer matrix multiplication
			startTime = System.nanoTime();
			matrix.divideAndConquer();
			endTime = System.nanoTime();
			// Print computation times in milliseconds
			System.out.println("Divide and Conquer:   " + (endTime - startTime));
	
			// Compute using Strassen algorithm for matrix multiplication
			startTime = System.nanoTime();
			matrix.strassenMultiply();
			endTime = System.nanoTime();
			// Print computation times in milliseconds
			System.out.println("Strassen Multiply:    " + (endTime - startTime));
			
			size *= 2;
			round++;
		}
	}

}
