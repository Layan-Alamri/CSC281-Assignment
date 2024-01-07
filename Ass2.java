package ass2;

//Assignment2 / CSC281

import java.util.*;
import java.math.*;

public class Ass2 {
    
    //--------------------------------------- Q1 ax=b(mod m). ---------------------------------------
    public static List<Integer> solveLinearCongruence(int a, int b, int m) {
        List<Integer> solutions = new ArrayList<>();

        int gcd = extendedEuclideanAlgorithm(a, m);
        if (b % gcd != 0) { 
            // No solutions exist
            return solutions;
        }

        int x0 = modularInverse(a / gcd, m / gcd);
        int x = (x0 * (b / gcd)) % (m / gcd);

        for (int i = 0; i < gcd; i++) {
            solutions.add(x + i * (m / gcd));
        }

        return solutions;
    }
    
    // Euclidean algorithm for finding the greatest common divisor (GCD)
    private static int extendedEuclideanAlgorithm(int a, int b) {
        if (b == 0) {
            return a;
        }
        return extendedEuclideanAlgorithm(b, a % b);
    }
    
    // Modular inverse using the extended Euclidean algorithm
    private static int modularInverse(int a, int m) {
        a = a % m;
        for (int x = 1; x < m; x++) {
            if ((a * x) % m == 1) {
                return x;
            }
        }
        return -1; // No modular inverse exists
    }
//--------------------------------------- Q2 Chinese Remainder ---------------------------------------
    
    public static long crt(int[] n, int[] l) {
        int len = n.length;

        // Compute the product of all moduli
        long prod = 1;
        for (int i = 0; i < len; i++) {
            prod *= n[i];
        }

        // Compute the array of factors
        long[] factors = new long[len];
        for (int i = 0; i < len; i++) {
            factors[i] = prod / n[i];
        }

        // Compute the array of inverses
        long[] inverses = new long[len];
        for (int i = 0; i < len; i++) {
            inverses[i] = modularInverse(factors[i], n[i]);
        }

        // Compute the solution using the Chinese Remainder Theorem
        long result = 0;
        for (int i = 0; i < len; i++) {
            result += (l[i] * factors[i] * inverses[i]) % prod;
        }
        result %= prod;
        return result;
    }

    //  method to compute the modular inverse of a modulo m
    private static long modularInverse(long a, long m) {
        long[] result = extendedEuclidean(a, m);
        long gcd = result[0];
        long x = result[1];
        if (gcd != 1) {
            throw new ArithmeticException("Modular inverse does not exist.");
        }
        return (x % m + m) % m;
    }

    //  method to compute the extended Euclidean algorithm
    private static long[] extendedEuclidean(long a, long b) {
        if (b == 0) {
            return new long[]{a, 1, 0};
        }
        long[] result = extendedEuclidean(b, a % b);
        long gcd = result[0];
        long x = result[2];
        long y = result[1] - (a / b) * result[2];
        return new long[]{gcd, x, y};
    }
    
    //--------------------------------------- Q3 Fermat Primality Test ---------------------------------------
    
     public static boolean fermat_Primality_Test(BigInteger big, int k) {
        if (big.equals(BigInteger.ONE)) {
            return false;  // 1 is not prime
        }

        Random random = new Random();

        for (int i = 0; i < k; i++) {
            BigInteger a = generateRandomBigInteger(BigInteger.TWO, big.subtract(BigInteger.ONE));
            BigInteger result = a.modPow(big.subtract(BigInteger.ONE), big);

            if (!result.equals(BigInteger.ONE)) {
                return false;  // n is composite
            }
        }

        return true;  // n is probably prime
    }

    //  method to generate a random BigInteger in the range [min, max]
    private static BigInteger generateRandomBigInteger(BigInteger min, BigInteger max) {
        BigInteger range = max.subtract(min).add(BigInteger.ONE);
        Random random = new Random();
        int bitLength = max.bitLength();

        BigInteger result;
        do {
            result = new BigInteger(bitLength, random);
        } while (result.compareTo(range) >= 0 || result.compareTo(BigInteger.ZERO) <= 0);

        return result.add(min);
    }

    //--------------------------------------- Q4 Modular Pow ---------------------------------------
    
   public static BigInteger modular_Pow(BigInteger base, BigInteger exponent, BigInteger modulus) {
        BigInteger result = BigInteger.ONE;
        base = base.mod(modulus);

        while (exponent.compareTo(BigInteger.ZERO) > 0) {
            if (exponent.testBit(0)) {
                result = result.multiply(base).mod(modulus);
            }
            base = base.multiply(base).mod(modulus);
            exponent = exponent.shiftRight(1);
        }

        return result;
    }

 //--------------------------------------- main ---------------------------------------
    
    public static void main(String[] args) {
        
    System.out.println(" \nThe solution for Q1: ");
        
        int a = 3;
        int b = 4;
        int m = 7;
        List<Integer> solutions = solveLinearCongruence(a, b, m);
       
        if (solutions.isEmpty()) {
            System.out.println("No solutions found.");
        }
        else {
            
            System.out.println(" Solutions:");
            for(int solution : solutions ){
                System.out.println(" " + a + "x modular congruence " + b + " (mod " + m + ") = " + solutions);
            }
        }
          
    System.out.println(" \nThe solution for Q2: ");
      
        int[] n = {3, 5, 7};  // List of moduli
        int[] l = {2, 3, 2};  // List of remainders

        long result = crt(n, l);
        System.out.println("Solution to the simultaneous congruences for x " + Arrays.toString(l) +
                " (mod " + Arrays.toString(n) + ") is =  " + result);
    
        
    System.out.println(" \nThe solution for Q3: ");
    
        BigInteger big = new BigInteger("7919");  // Number to be tested
             int k = 5;  // Number of iterations

             boolean isPrime = fermat_Primality_Test(big, k);
             System.out.println(big + " is " + (isPrime ? "probably prime." : "composite."));

    System.out.println(" \nThe solution for Q4: ");
         
         BigInteger base = new BigInteger("3");
             BigInteger exponent = new BigInteger("644");
             BigInteger modulus = new BigInteger("645");

             BigInteger result1 = modular_Pow(base, exponent, modulus);
             BigInteger result2 = base.modPow(exponent, modulus);

              System.out.println("(base^exponent) mod modulus = " + result1);
              System.out.println("base.modPow(exponent, modulus) = " + result2);
              System.out.println("Results are equal: " + result1.equals(result2));
              
    }// End main 
}// End class


 