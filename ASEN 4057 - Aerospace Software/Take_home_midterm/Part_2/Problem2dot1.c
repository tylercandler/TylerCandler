#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int BinarySearch(int arr[], int elements, int target){
  //Inputs: Array, number of indices, and target number
  //Output: index of targeted number in the associated array
  int right = elements - 1; //right side index
  int left = 0; //left side index

  while(right >= left){
    int middle = (right + left)/2; // find middle index
    if (arr[middle] == target){ // if middle element of array is target
      return middle; // return middle index
    }
    if (arr[middle] > target){ // if middle element is greater than target
      right = middle - 1; //only continue search on left side of array
    }
    if (arr[middle] < target){ //if middle element is less than target
      left = middle + 1; //only continue search on right side of array
    }
  }
  return -50; // if target does not exist return arbitray negative number
}

int LinearSearch(int array[], int elements, int target){
  //Inputs: Array, number of indices, and target number
  //Output: index of targeted number in the associated array
  for (int i = 0; i <= (elements-1); i+=1){ //for loop through all indices of array
    if (array[i] == target){ //stop once target is found
        return i;
    }
  }
  return -50;
}


int main( int argc, char *argv[] ){
  int array[100]; //declare array as an overestimate
  int index; //declare index
  int target; //declare target
  //printf("File number is %c\n", argv[1][5]);
  char file_char = argv[1][5];
  int file_num;
  sscanf(&argv[1][5],"%d",&file_num);// use sscanf to convert argv[2] as an integer to target
  sscanf(argv[2],"%d",&target); // same thing with file number, except double extracting from char*
  int i = 0; //initialize counter
  FILE *infile; //file pointer declaration
  int elements = 0; // initialize elements
  printf("Examining %s\n", argv[1]);
  if (infile = (fopen(argv[1], "r"))){ // if file exists, open it
    while (fscanf(infile, "%d", &array[i]) != EOF){ //scan through each line of the file
      i++; // increase counter
    }
    fclose(infile); // close file
    array[i] = '\0'; //declare array of numbers, for now terminate with null

    for (i = 0; array[i] != '\0'; i++){
      // printf("%d\n", array[i]);
      elements+=1;
    }
  } else {printf("The file %s does not exist\n", argv[1]);} //print if file does not exist

if (file_num % 2 == 0){ // if file number is even (sorted)
  index = BinarySearch(array, elements, target);//do BinarySearch
}
else{//if file number is odd (unsorted)
  index = LinearSearch(array, elements, target); // do LinearSearch
}


char output_filename[]  = "ArrayP.out";
output_filename[5] = file_char;



//printf("OUTPUT %s\n", output_filename);



//
FILE *outfile = fopen(output_filename,"ab+"); // open output file and create one if it doesn't already exist
if(index != -50){//if target does not exist in array
  fprintf(outfile, "The number %d is at index %d in %s\n", target, index, argv[1]);
  // printf("The number %d is at index %d in %s\n", target, index, argv[1]);
}else{
  fprintf(outfile, "The number %d does not exist in %s\n", target, argv[1]);
  // printf("The number %d does not exist in %s\n", target, argv[1]);
}
fclose(outfile);


return 0;
}







// -----------------------------------------------------------------------------
