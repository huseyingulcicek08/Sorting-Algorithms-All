public class TestSortingAlgorithms {
/*
    DO NOT WRITE ANY CODE HERE

    -------SAMPLE OUTPUT-1 --------------

    Bubble Sort - TRUE
    Selection Sort - TRUE
    Insertion Sort - TRUE
    Merge Sort - TRUE
    Quick Sort - TRUE
    Heap Sort - TRUE
    Counting Sort - TRUE
    Radix Sort - TRUE

    *****Summary of Sorting Algorithms*****
    Using Bubble Sort It took->17148 milliseconds
    Using Selection Sort It took->4621 milliseconds
    Using Insertion Sort It took->941 milliseconds
    Using Merge Sort It took->16 milliseconds
    Using Quick Sort It took->31 milliseconds
    Using Heap Sort It took->16 milliseconds
    Using Counting Sort It took->15 milliseconds
    Using Radix Sort It took->32 milliseconds

    Acording to the given imput the most successful algorithm is Counting Sort

    Your Grade is->90


	-----------------------------------------------------------------------------------


    -------------SAMPLE OUTPUT-2-------------

    Bubble Sort - TRUE
    Selection Sort - TRUE
    Insertion Sort - TRUE
    Merge Sort - TRUE
    Quick Sort - TRUE
    Heap Sort - TRUE
    Counting Sort - TRUE
    Radix Sort - TRUE

    *****Summary of Sorting Algorithms*****
    Using Bubble Sort It took->16148 milliseconds
    Using Selection Sort It took->5121 milliseconds
    Using Insertion Sort It took->891 milliseconds
    Using Merge Sort It took->21 milliseconds
    Using Quick Sort It took->27 milliseconds
    Using Heap Sort It took->18 milliseconds
    Using Counting Sort It took->15 milliseconds
    Using Radix Sort It took->15 milliseconds

    Acording to the given imput the most successful algorithms are Counting Sort and Radix Sort

    Your Grade is->90

	-----------------------------------------------------------------------------------

    NOTE: WRITING 8 ALGORITHMS CORRECTLY WORTH 90 POINTS,
    REMAINING 10 POINTS WILL BE GIVEN BY TEACHING ASISTANT ACORDING TO YOUR
    public void summaryofSortingAlgorithms() METHOD.

    */



    public static void main(String[] args) {
        int grade=0;
        int A[]=new int[100000];

        SortingAlgortithms sort=new SortingAlgortithms();
        shuffleArray(A);
        sort.BubbleSort(A);

        if (isSorted(A)){
            System.out.println("Bubble Sort - TRUE");
            grade += 5;

        }
        else{
            System.out.println("Bubble Sort - FALSE");
        }
        shuffleArray(A);
        sort.SelectionSort(A);
        if (isSorted(A)){
            System.out.println("Selection Sort - TRUE");
            grade += 5;
        }
        else{
            System.out.println("Selection Sort - FALSE");
        }
        shuffleArray(A);
        sort.InsertionSort(A);
        if (isSorted(A)){
            System.out.println("Insertion Sort - TRUE");
            grade += 7;
        }
        else{
            System.out.println("Insertion Sort - FALSE");
        }
        shuffleArray(A);
        sort.MergeSort(A);
        if (isSorted(A)){
            System.out.println("Merge Sort - TRUE");
            grade += 20;
        }
        else{
            System.out.println("Merge Sort - FALSE");
        }
        shuffleArray(A);
        sort.QuickSort(A);
        if (isSorted(A)){
            System.out.println("Quick Sort - TRUE");
            grade += 20;
        }
        else{
            System.out.println("Quick Sort - FALSE");
        }
        shuffleArray(A);
        sort.HeapSort(A);
        if (isSorted(A)){
            System.out.println("Heap Sort - TRUE");
            grade += 14;
        }
        else{
            System.out.println("Heap Sort - FALSE");
        }
        shuffleArray(A);
        sort.CountingSort(A);

        if (isSorted(A)){
            System.out.println("Counting Sort - TRUE");
            grade += 7;
        }
        else{
            System.out.println("Counting Sort - FALSE");
        }
        shuffleArray(A);
        sort.RadixSort(A);
        if (isSorted(A)){
            System.out.println("Radix Sort - TRUE");
            grade += 12;
        }
        else{
            System.out.println("Radix Sort - FALSE");
        }
        System.out.println();
        sort.summaryofSortingAlgorithms();
        System.out.println();
        System.out.println("Your Grade is->"+grade);
    }

    public static void shuffleArray(int [] a){
        for(int i=0;i<a.length;i++)
            a[i]=(int)(Math.random()*10000);
    }

    public static boolean isSorted(int a[]){
        boolean sorted=true;
        for(int i=0;i<a.length-1;i++)
            if(a[i]>a[i+1])
                sorted=false;
        return sorted;
    }
}
