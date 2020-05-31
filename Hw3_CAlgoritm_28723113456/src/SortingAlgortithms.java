import java.util.Arrays;
import java.io.*;
import java.util.*;
import java.util.concurrent.TimeUnit;

public class SortingAlgortithms {
    private static SortingAlgortithms RadixSort;
    private int[] temp_array;

/*  implement sorting algorithms below
•	Bubble Sort
•	Selection Sort
•	Insertion Sort
•	Merge Sort
•	Quick Sort
•	Heap Sort
•	Counting Sort
•	Radix Sort

   Fill in the method.
   Implement sorting algorithms.
   You can use extra method.

*/

    public void BubbleSort(int A[]) {
        int n = A.length;
        for (int i = 0; i < n-1; i++)
            for (int j = 0; j < n-i-1; j++)
                if (A[j] > A[j+1])
                {
                    // swap arr[j+1] and arr[i]
                    int temp = A[j];
                    A[j] = A[j+1];
                    A[j+1] = temp;
                }
    }

    public void SelectionSort(int A[]) {
        int n = A.length;

        // One by one move boundary of unsorted subarray
        for (int i = 0; i < n-1; i++)
        {
            // Find the minimum element in unsorted array
            int min_idx = i;
            for (int j = i+1; j < n; j++)
                if (A[j] < A[min_idx])
                    min_idx = j;

            // Swap the found minimum element with the first
            // element
            int temp = A[min_idx];
            A[min_idx] = A[i];
            A[i] = temp;
        }

    }

    public void InsertionSort(int A[]) {
        int n = A.length;
        for (int i = 1; i < n; ++i) {
            int key = A[i];
            int j = i - 1;

            /* Move elements of arr[0..i-1], that are
               greater than key, to one position ahead
               of their current position */
            while (j >= 0 && A[j] > key) {
                A[j + 1] = A[j];
                j = j - 1;
            }
            A[j + 1] = key;
        }

    }

    public void MergeSort(int A[]) {
        if(A == null)
        {
            return;
        }

        if(A.length > 1)
        {
            int mid = A.length / 2;

            // Split left part
            int[] left = new int[mid];
            for(int i = 0; i < mid; i++)
            {
                left[i] = A[i];
            }

            // Split right part
            int[] right = new int[A.length - mid];
            for(int i = mid; i < A.length; i++)
            {
                right[i - mid] = A[i];
            }
            MergeSort(left);
            MergeSort(right);

            int i = 0;
            int j = 0;
            int k = 0;

            // Merge left and right arrays
            while(i < left.length && j < right.length)
            {
                if(left[i] < right[j])
                {
                    A[k] = left[i];
                    i++;
                }
                else
                {
                    A[k] = right[j];
                    j++;
                }
                k++;
            }
            // Collect remaining elements
            while(i < left.length)
            {
                A[k] = left[i];
                i++;
                k++;
            }
            while(j < right.length)
            {
                A[k] = right[j];
                j++;
                k++;
            }
        }

    }

    public void QuickSort(int A[]) {
        if (A == null || A.length == 0) {
            return;
        }
        this.temp_array = A;
        int len = A.length;
        quickSort(0, len - 1);
    }
    private void quickSort(int low_index, int high_index) {

        int i = low_index;
        int j = high_index;
        // calculate pivot number
        int pivot = temp_array[low_index+(high_index-low_index)/2];
        // Divide into two arrays
        while (i <= j) {
            while (temp_array[i] < pivot) {
                i++;
            }
            while (temp_array[j] > pivot) {
                j--;
            }
            if (i <= j) {
                exchangeNumbers(i, j);
                //move index to next position on both sides
                i++;
                j--;
            }
        }
        // call quickSort() method recursively
        if (low_index < j)
            quickSort(low_index, j);
        if (i < high_index)
            quickSort(i, high_index);
    }

    private void exchangeNumbers(int i, int j) {
        int temp = temp_array[i];
        temp_array[i] = temp_array[j];
        temp_array[j] = temp;
    }

    public void HeapSort(int A[]){
        int N = A.length;

        // Build heap (rearrange array)
        for (int i = N / 2 - 1; i >= 0; i--) {
            heapify(A, N, i);
        }
        // One by one extract an element from heap
        for (int i = N - 1; i >= 0; i--) {
            // Move current root to end
            int temp = A[0];
            A[0] = A[i];
            A[i] = temp;

            // Call max heapify on the reduced heap
            heapify(A, i, 0);
        }

    }
    public void heapify(int[] numbers, int n, int i) {
        int largest = i; // initialize largest as root
        int l = 2*i + 1; // left = 2*i + 1
        int r = 2*i + 2; // right = 2*i + 2

        // if left child is larger than root
        if (l < n && numbers[l] > numbers[largest]) {
            largest = l;
        }

        // if right child is larger than largest so far
        if (r < n && numbers[r] > numbers[largest]) {
            largest = r;
        }

        // if largest is not root
        if (largest != i) {
            int temp = numbers[i];
            numbers[i] = numbers[largest];
            numbers[largest] = temp;

            // Recursively heapify the affected sub-tree
            heapify(numbers, n, largest);
        }
    }

    public void CountingSort( int A[]){
        int max = Arrays.stream(A).max().getAsInt();
        int min = Arrays.stream(A).min().getAsInt();
        int range = max - min + 1;
        int count[] = new int[range];
        int output[] = new int[A.length];
        for (int i = 0; i < A.length; i++)
        {
            count[A[i] - min]++;
        }

        for (int i = 1; i < count.length; i++)
        {
            count[i] += count[i - 1];
        }

        for (int i = A.length - 1; i >= 0; i--)
        {
            output[count[A[i] - min] - 1] = A[i];
            count[A[i] - min]--;
        }

        for (int i = 0; i < A.length; i++)
        {
            A[i] = output[i];
        }
    }




        public void RadixSort(int[] A)
        {
            RadixSort.sort(A, 10);
        }
    public static void sort(int[] array, int radix) {
        if (array.length == 0) {
            return;
        }

        // Determine minimum and maximum values
        int minValue = array[0];
        int maxValue = array[0];
        for (int i = 1; i < array.length; i++) {
            if (array[i] < minValue) {
                minValue = array[i];
            } else if (array[i] > maxValue) {
                maxValue = array[i];
            }
        }

        // Perform counting sort on each exponent/digit, starting at the least
        // significant digit
        int exponent = 1;
        while ((maxValue - minValue) / exponent >= 1) {
            RadixSort.countingSortByDigit(array, radix, exponent, minValue);
            exponent *= radix;
        }
    }

    private static void countingSortByDigit(
            int[] array, int radix, int exponent, int minValue) {
        int bucketIndex;
        int[] buckets = new int[radix];
        int[] output = new int[array.length];

        // Initialize bucket
        for (int i = 0; i < radix; i++) {
            buckets[i] = 0;
        }

        // Count frequencies
        for (int i = 0; i < array.length; i++) {
            bucketIndex = (int)(((array[i] - minValue) / exponent) % radix);
            buckets[bucketIndex]++;
        }

        // Compute cumulates
        for (int i = 1; i < radix; i++) {
            buckets[i] += buckets[i - 1];
        }

        // Move records
        for (int i = array.length - 1; i >= 0; i--) {
            bucketIndex = (int)(((array[i] - minValue) / exponent) % radix);
            output[--buckets[bucketIndex]] = array[i];
        }

        // Copy back
        for (int i = 0; i < array.length; i++) {
            array[i] = output[i];
        }
    }

    public void summaryofSortingAlgorithms(){
        System.out.println("*****Summary of Sorting Algorithms*****");
//WRITE YOUR SUMMARY HERE

                int[] anArray;
                anArray = new int[1000];
                Random generator = new Random();
                for(int i=0; i<1000; i++){
                    anArray[i] = (generator.nextInt(1000)+1);
                }
                Date before = new Date();
                Date after;
                Arrays.sort(anArray);
                after = new Date();
                System.out.println(after.getTime()-before.getTime());

                System.out.println(Arrays.toString(anArray));
        long startingTime=System.nanoTime();
        Arrays.sort(anArray);
        long endTime=System.nanoTime();
        System.out.println("Sorting time: "+(endTime-startingTime)+"ns");


    }
}
