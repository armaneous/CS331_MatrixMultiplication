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

import java.util.*;

public class Matrix {
	private int[][] resultMatrix;			// Resultant matrix after multiplication
	private int[][] a;						// Initial matrix created
	private int[][] b;						// Second matrix created
	private Random rand = new Random();		// For randomly generating numbers
	
	/**
	 * Constructor for Matrix, initializes two matrices
	 * and instantiates size of resultant matrix
	 * @param size		>>	Size of matrices to be generated
	 */
	public Matrix(int size){
		a = generateRandomMatrix(size);
		b = generateRandomMatrix(size);
		resultMatrix = new int[size][size];
	}
	
	/**
	 * Generate random numbers from 0 to 10, exclusive,
	 * and fill arrays with randomly generated numbers
	 * @param size		>>	Size of matrices to be generated
	 * @return			>>	Returns a matrix to be used
	 */
	private int[][] generateRandomMatrix(int size){
		int[][] randMatrix = new int[size][size];
		for(int i = 0; i < size; i++)
			for(int j = 0; j < size; j++)
				randMatrix[i][j] = rand.nextInt(10);
		return randMatrix;
	}
	
	/**
	 * Classical matrix multiplication, which runs to the
	 * order of O(n^3) time
	 * @return		>>	Returns resultant matrix
	 */
	public int[][] naiveMultiply(){
		for(int i = 0; i < resultMatrix.length; i++)
			for(int j = 0; j < resultMatrix.length; j++){
				resultMatrix[i][j] = 0;
				for(int k = 0; k < resultMatrix.length; k++)
					resultMatrix[i][j] += a[i][k] * b[k][j];
			}
		return resultMatrix;
	}
	
	/**
	 * Divide and Conquer method for multiplying matrices, uses
	 * private helper method
	 * @return		>>	Returns resultant matrix
	 */
	public int[][] divideAndConquer(){
		return divideAndConquer(a, b);
	}
	
	/**
	 * Private helper method for Divide and Conquer method of
	 * matrix multiplication
	 * @param matrixA	>>	First matrix input for multiplication
	 * @param matrixB	>>	Second matrix input for multiplication
	 * @return			>>	Returns merged matrix 
	 */
	private int[][] divideAndConquer(int[][] matrixA, int[][] matrixB){
        int[][] c11, c12, c21, c22;
	    if (matrixA.length == 2){
	         // Calculate and return base case
	    	int[][] c = new int[2][2];
	    	c[0][0] = (matrixA[0][0] * matrixB[0][0]) + (matrixA[0][1] * matrixB[1][0]);
	    	c[0][1] = (matrixA[0][0] * matrixB[0][1]) + (matrixA[0][1] * matrixB[1][1]);
	    	c[1][0] = (matrixA[1][0] * matrixB[0][0]) + (matrixA[1][1] * matrixB[1][0]);
	    	c[1][1] = (matrixA[1][0] * matrixB[0][1]) + (matrixA[1][1] * matrixB[1][1]);
	    	return c;
	    }
	    else {
			// Partition matrices into four equal parts
	    	int[][][] mArrayA = partitionMatrix(matrixA);
	    	int[][][] mArrayB = partitionMatrix(matrixB);
	    	// Each fourth repeats to divide and conquer
	        c11 = addMatrix(divideAndConquer(mArrayA[0],mArrayB[0]),divideAndConquer(mArrayA[1],mArrayB[2]));
	        c12 = addMatrix(divideAndConquer(mArrayA[0],mArrayB[1]),divideAndConquer(mArrayA[1],mArrayB[3]));
	        c21 = addMatrix(divideAndConquer(mArrayA[2],mArrayB[0]),divideAndConquer(mArrayA[3],mArrayB[2]));
	        c22 = addMatrix(divideAndConquer(mArrayA[2],mArrayB[1]),divideAndConquer(mArrayA[3],mArrayB[3]));
	        // Merge matrices into one, full, matrix
	        return mergeMatrix(c11, c12, c21, c22);
	    }
	}
	
	/**
	 * Uses the Strassen method for matrix multiplication,
	 * which requires 7 multiplications and 18 additions
	 * or subtractions. Utilizes private helper method.
	 * @return		>>	Returns resultant matrix
	 */
	public int[][] strassenMultiply(){
		return strassenMultiply(a.length, a, b, resultMatrix);
	}
	
	/**
	 * Private helper method for Strassen method
	 * @param size		>>	Takes in size of arrays fed into method
	 * @param matrixA	>>	First matrix used for multiplication
	 * @param matrixB	>>	Second matrix used for multiplication
	 * @param result	>>	Matrix for which the output should be written to
	 * @return			>>	Returns merged matrix
	 */
	private int[][] strassenMultiply(int size, int[][] matrixA, int[][] matrixB, int[][] result){
		int[][] p, q, r, s, t, u, v;
		int[][] c11, c12, c21, c22;
	    if (matrixA.length == 2){
	         // Calculate and return base case
	    	result[0][0] = (matrixA[0][0] * matrixB[0][0]) + (matrixA[0][1] * matrixB[1][0]);
	    	result[0][1] = (matrixA[0][0] * matrixB[0][1]) + (matrixA[0][1] * matrixB[1][1]);
	    	result[1][0] = (matrixA[1][0] * matrixB[0][0]) + (matrixA[1][1] * matrixB[1][0]);
	    	result[1][1] = (matrixA[1][0] * matrixB[0][1]) + (matrixA[1][1] * matrixB[1][1]);
	    	return result;
	    }
	    else {
			// Partition matrices into four equal parts
	    	int[][][] mArrayA = partitionMatrix(matrixA);
	    	int[][][] mArrayB = partitionMatrix(matrixB);
	    	
	    	p = q = r = s = t = u = v = new int[mArrayA[0].length][mArrayA[0].length];
	        strassenMultiply(size/2, addMatrix(mArrayA[0],mArrayA[3]), addMatrix(mArrayB[0],mArrayB[3]), p);
	        strassenMultiply(size/2, addMatrix(mArrayA[2],mArrayA[3]), mArrayB[0], q);
	        strassenMultiply(size/2, mArrayA[0], subMatrix(mArrayB[1],mArrayB[3]), r);
	        strassenMultiply(size/2, mArrayA[3], subMatrix(mArrayB[2],mArrayB[1]), s);
	        strassenMultiply(size/2, addMatrix(mArrayA[0],mArrayA[1]), mArrayB[3], t);
	        strassenMultiply(size/2, subMatrix(mArrayA[2],mArrayA[0]), addMatrix(mArrayB[0],mArrayB[1]), u);
	        strassenMultiply(size/2, subMatrix(mArrayA[1],mArrayA[3]), addMatrix(mArrayB[2],mArrayB[3]), v);
	        c11 = subMatrix(addMatrix(p, s), addMatrix(t, v));
	        c12 = addMatrix(r, t);
	        c21 = addMatrix(q, s);
	        c22 = subMatrix(addMatrix(p, r), addMatrix(q, u));
	        // Merge matrices into one, full, matrix
	        return mergeMatrix(c11, c12, c21, c22);
	    }
	}
	
	/**
	 * Method used for adding two matrices of equal size
	 * @param matrixA	>>	First matrix used for adding
	 * @param matrixB	>>	Second matrix used for adding
	 * @return			>>	Returns sum of two matrices
	 */
	private int[][] addMatrix(int[][] matrixA, int[][] matrixB) {
		int[][] result = new int[matrixA.length][matrixA.length];
		for(int i = 0; i < matrixA.length; i++)
			for(int j = 0; j < matrixA.length; j++) {
				result[i][j] = matrixA[i][j] + matrixB[i][j];
			}
		return result;
	}

	/**
	 * Method used for subtracting two matrices of equal size
	 * @param matrixA	>>	First matrix used for subtracting
	 * @param matrixB	>>	Second matrix used for subtracting
	 * @return			>>	Returns difference of two matrices
	 */
	private int[][] subMatrix(int[][] matrixA, int[][] matrixB) {
		int[][] result = new int[matrixA.length][matrixA.length];
		for(int i = 0; i < matrixA.length; i++)
			for(int j = 0; j < matrixA.length; j++) {
				result[i][j] = matrixA[i][j] - matrixB[i][j];
			}
		return result;
	}

	/**
	 * Partitions a square matrix of size that is to the power of 2
	 * into four equal sub-matrices of original matrix
	 * @param matrix	>>	Initial matrix to be partitioned
	 * @return			>>	Returns array of size 4 with resulting quadrants
	 */
	private int[][][] partitionMatrix(int[][] matrix){
		int[][] matrix11 = new int[matrix.length/2][matrix.length/2];
		int[][] matrix12 = new int[matrix.length/2][matrix.length/2];
		int[][] matrix21 = new int[matrix.length/2][matrix.length/2];
		int[][] matrix22 = new int[matrix.length/2][matrix.length/2];
		for(int i = 0; i < matrix.length/2; i++)
			for(int j = 0; j < matrix.length/2; j++){
				matrix11[i][j] = matrix[i][j];
				matrix12[i][j] = matrix[i][j + matrix.length/2];
				matrix21[i][j] = matrix[i + matrix.length/2][j];
				matrix22[i][j] = matrix[i + matrix.length/2][j + matrix.length/2];
			}
		return new int[][][] { matrix11, matrix12, matrix21, matrix22};
	}
	
	/**
	 * Merges four quadrants, or sub-matrices, into one complete matrix
	 * @param matrix11	>>	First quadrant, or top left, of full matrix
	 * @param matrix12	>>	Second quadrant, or top right, of full matrix
	 * @param matrix21	>>	Third quadrant, or bottom left, of full matrix
	 * @param matrix22	>>	Fourth quadrant, or bottom right, of full matrix
	 * @return			>>	Returns final, merged, matrix
	 */
	private int[][] mergeMatrix(int[][] matrix11, int[][] matrix12, int[][] matrix21, int[][] matrix22){
		int[][] result = new int[matrix11.length*2][matrix11.length*2];
		for(int i = 0; i < matrix11.length; i++)
			for(int j = 0; j < matrix11.length; j++){
				result[i][j] = matrix11[i][j];
				result[i][j + matrix12.length] = matrix12[i][j];
				result[i + matrix21.length][j] = matrix21[i][j];
				result[i + matrix22.length][j + matrix22.length] = matrix22[i][j];
			}
		return result;
	}
	
	/**
	 * Converts resultant matrix into a neatly formatted string to
	 * when printed
	 */
	public String toString(){
		String s = "";
		for(int i = 0; i < resultMatrix.length; i++)
			for(int j = 0; j < resultMatrix.length; j++){
				s += String.valueOf(resultMatrix[i][j]) + " \t";
				if(j == resultMatrix.length - 1)
					s += "\n";
			}
		return s;
	}
}
