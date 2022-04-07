/*
 Copyright (c) 2019, Shinichi Morishita
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 The views and conclusions contained in the software and documentation are those
 of the authors and should not be interpreted as representing official policies,
 either expressed or implied, of the FreeBSD Project.
 */

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mTR.h"
#include <time.h>


void free_global_variables_and_exit(){
    // If any of the above global variables failed to be allocated in the heap, free other variables and exit.
    if(nextReadID        != NULL){ free(nextReadID); }
    if(pow4              != NULL){ free(pow4); }
    if(orgInputString    != NULL){ free(orgInputString); }
    if(inputString       != NULL){ free(inputString); }
    if(inputString_w_rand!= NULL){ free(inputString_w_rand); }
    if(count             != NULL){ free(count); }
    if(freqNode          != NULL){ free(freqNode); }
    for(int i = 0; i < PrimeMax; i++)
        if( gaps[i] != NULL ){ free(freqNode[i]); }
    if(sortedString      != NULL){ free(sortedString); }
    if(directional_index_tmp != NULL){ free(directional_index_tmp); }
    if(directional_index != NULL){ free(directional_index); }
    if(directional_index_end != NULL){ free(directional_index_end); }
    if(directional_index_w != NULL){ free(directional_index_w); }
    if(vector0           != NULL){ free(vector0); }
    if(vector1           != NULL){ free(vector1); }
    if(vector2           != NULL){ free(vector2); }
    if(freq_interval_len != NULL){ free(freq_interval_len); }
    if(WrapDP            != NULL){ free(WrapDP); }
    if(alignment_input   != NULL){ free(alignment_input); }
    if(alignment_symbols != NULL){ free(alignment_symbols); }
    if(alignment_repeats != NULL){ free(alignment_repeats); }
    for(int i=0; i<(MAX_PERIOD + 1); i++)
        if( consensus[i] != NULL ){ free(consensus[i]); }
    if(consensus         != NULL){ free(consensus); }
    for(int i=0; i<(MAX_PERIOD + 1); i++)
        if( gaps[i] != NULL ){ free(gaps[i]); }
    if( gaps             != NULL ){ free(gaps); }
    fprintf(stderr, "cannot allocate space for one of global variables in the heap.\n");
    exit(EXIT_FAILURE);
}

void malloc_global_variables(){
    
    // Allocate the main memory for global variables in the heap
    nextReadID = (char *)malloc(sizeof(char) * MAX_ID_LENGTH);
    if( nextReadID == NULL ){ free_global_variables_and_exit(); }
    pow4  = (int *)malloc(sizeof(int) * (maxKmer+1));
    if( pow4 == NULL ){ free_global_variables_and_exit(); }
    pow4[0] = 1;
    for(int i=1; i<= maxKmer; i++)
        pow4[i] = pow4[i-1] * 4;
    orgInputString  = (int *)malloc(sizeof(int) * MAX_INPUT_LENGTH);
    if( orgInputString == NULL ){ free_global_variables_and_exit(); }
    inputString     = (int *)malloc(sizeof(int) * MAX_INPUT_LENGTH);
    if( inputString == NULL ){ free_global_variables_and_exit(); }
    inputString_w_rand = (int *)malloc(sizeof(int) * MAX_INPUT_LENGTH);
    if( inputString_w_rand == NULL ){ free_global_variables_and_exit(); }
    count           = (int *)malloc( sizeof(int) * pow4[count_maxKmer]);
    if( count == NULL ){ free_global_variables_and_exit(); }
    freqNode = malloc(sizeof(int *) * PrimeMax);
    if( freqNode == NULL ){ free_global_variables_and_exit(); }
    for(int i = 0; i < PrimeMax; i++){
        freqNode[i] = malloc(sizeof(int) * 2);
        if( freqNode[i] == NULL ){ free_global_variables_and_exit(); }
    }
    sortedString    = (int *)malloc(sizeof(int) * MAX_INPUT_LENGTH);
    if( sortedString == NULL ){ free_global_variables_and_exit(); }
    directional_index_tmp   = (double *)malloc(sizeof(double) * MAX_INPUT_LENGTH);
    if( directional_index_tmp == NULL ){ free_global_variables_and_exit(); }
    directional_index       = (double *)malloc(sizeof(double) * MAX_INPUT_LENGTH);
    if( directional_index == NULL ){ free_global_variables_and_exit(); }
    directional_index_end   = (int *)malloc(sizeof(int) * MAX_INPUT_LENGTH);
    if( directional_index_end == NULL ){ free_global_variables_and_exit(); }
    directional_index_w     = (int *)malloc(sizeof(int) * MAX_INPUT_LENGTH);
    if( directional_index_w == NULL ){ free_global_variables_and_exit(); }
    vector0 = (int *)malloc(sizeof(int) * 4 * BLK);
    if( vector0 == NULL ){ free_global_variables_and_exit(); }
    vector1 = (int *)malloc(sizeof(int) * 4 * BLK);
    if( vector1 == NULL ){ free_global_variables_and_exit(); }
    vector2 = (int *)malloc(sizeof(int) * 4 * BLK);
    if( vector2 == NULL ){ free_global_variables_and_exit(); }
    freq_interval_len   = (double *)malloc( sizeof(double) * MAX_INPUT_LENGTH);
    if( freq_interval_len == NULL ){ free_global_variables_and_exit(); }
    WrapDP          = (int *)malloc(sizeof(int) * WrapDPsize);
    if( WrapDP == NULL ){ free_global_variables_and_exit(); }
    alignment_input = (char *)malloc(sizeof(char) * MAX_INPUT_LENGTH);
    if( alignment_input == NULL ){ free_global_variables_and_exit(); }
    alignment_symbols = (char *)malloc(sizeof(char) * MAX_INPUT_LENGTH);
    if( alignment_symbols == NULL ){ free_global_variables_and_exit(); }
    alignment_repeats = (char *)malloc(sizeof(char) * MAX_INPUT_LENGTH);
    if( alignment_repeats == NULL ){ free_global_variables_and_exit(); }
    consensus = malloc(sizeof(int *) * (MAX_PERIOD + 1));
    if( consensus == NULL ){ free_global_variables_and_exit(); }
    for(int i=0; i<(MAX_PERIOD + 1); i++){
        consensus[i] = malloc(sizeof(int) * 5);
        if( consensus[i] == NULL ){ free_global_variables_and_exit(); }
    }
    gaps = malloc(sizeof(int *) * (MAX_PERIOD + 1));
    if( gaps == NULL ){ free_global_variables_and_exit(); }
    for(int i=0; i<(MAX_PERIOD + 1); i++){
        gaps[i] = malloc(sizeof(int) * 4);
        if( gaps[i] == NULL ){ free_global_variables_and_exit(); }
    }
#ifdef DEBUG
    fprintf(stderr, "Succeeded in malooc (handle_one_file.c)\n");
#endif
}

void free_global_variables(){
    free(nextReadID);
    free(pow4);
    free(orgInputString);
    free(inputString);
    free(inputString_w_rand);
    free(count);
    for(int i=0; i<PrimeMax; i++){ free(freqNode[i]); }
    free(freqNode);
    free(sortedString);
    free(directional_index_tmp);
    free(directional_index);
    free(directional_index_end);
    free(directional_index_w);
    free(vector0);
    free(vector1);
    free(vector2);
    free(freq_interval_len);
    free(WrapDP);
    free(alignment_input);
    free(alignment_symbols);
    free(alignment_repeats);
    for(int i=0; i<(MAX_PERIOD + 1); i++){ free(consensus[i]); }
    for(int i=0; i<(MAX_PERIOD + 1); i++){ free(gaps[i]); }
    free(consensus);
    free(gaps);
#ifdef DEBUG
    fprintf(stderr, "Succeeded in free (handle_one_file.c)\n");
#endif
}

int char2int(char c){
    int charCode;
    switch(c){
        case 'A':
        case 'a':
            charCode = 0; break;
        case 'C':
        case 'c':
            charCode = 1; break;
        case 'G':
        case 'g':
            charCode = 2; break;
        case 'T':
        case 't':
            charCode = 3; break;
        default:
            fprintf(stderr, "Invalid character: %c \n", c); exit(EXIT_FAILURE);
    }
    return(charCode);
}

FILE* init_handle_one_file(char *inputFile){
    FILE *fp = fopen(inputFile, "r");
    if(fp == NULL){
        fprintf(stderr, "fatal error: cannot open %s\n", inputFile);
        fflush(stderr);
        exit(EXIT_FAILURE);
    }
    tmp_read_cnt = -1;
    return(fp);
}

void return_one_read(FILE *fp, Read *currentRead){
    // Return 0 if it is the last read
    char *s = (char *)malloc(sizeof(char)*BLK);
    int i;
    int cnt=0;
    int no_read = 1;
    
    while (fgets(s, BLK, fp) != NULL) { // Feed a string of size BLK from fp into string s
        no_read = 0;
        
        if(s[0] == '>'){
            if(tmp_read_cnt == -1){ // This is the first read
                tmp_read_cnt = 0;
                // Feed the ID of the current read into the ID of nextRead
                for(i=1; s[i]!='\0' && s[i]!='\n' && s[i]!='\r' && i<BLK; i++)
                    nextReadID[i-1] = s[i];
                nextReadID[i-1] = '\0';
                // Move on to feed the current string
            }else{
                // Set the ID of currentRead to the ID of nextRead
                int j;
                for(j=0;  nextReadID[j] != '\0'; j++)
                    currentRead->ID[j] = nextReadID[j];
                currentRead->ID[j] = '\0';
                
                // Feed the ID of the current read into the ID of nextRead
                for(i=1; s[i]!='\0' && s[i]!='\n' && s[i]!='\r' && i<BLK; i++)
                    nextReadID[i-1] = s[i];
                nextReadID[i-1] = '\0';
                
                // Finalize the currentRead string by appending '\0'
                //currentRead.string[cnt] = '\0';
                currentRead->len = cnt;
                tmp_read_cnt++;
                free(s);
                return;
                //return(currentRead);
            }
        }else{
            // Feed the string
            for(i=0; s[i]!='\0' && s[i]!='\n' && s[i]!='\r'; i++){
            //for(i=0; s[i]!='\0' && s[i]!='\n'; i++){
                currentRead->codedString[cnt++] = char2int(s[i]);
                if( MAX_INPUT_LENGTH <= cnt){
                    fprintf(stderr, "fatal error: The length %d is tentatively at most %i.\nread ID = %s\nSet MAX_INPUT_LENGTH to a larger value", cnt, MAX_INPUT_LENGTH, currentRead->ID);
                    free_global_variables_and_exit();
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
    free(s);
    if(no_read == 1){  // No reads
        currentRead->len = 0;
        //return(currentRead);
    }else{
        // Process the last read.
        // Set the ID of currentRead to the ID of nextRead
        int j;
        for(j=0;  nextReadID[j] != '\0'; j++)
            currentRead->ID[j] = nextReadID[j];
        currentRead->ID[j] = '\0';
        // Finalize the currentRead string by appending '\0'
        //currentRead.string[cnt] = '\0';
        currentRead->len = cnt;
        tmp_read_cnt++;
        //return(currentRead); // The last read
    }
}

int handle_one_file(char *inputFile, int print_alignment, int numprocs, int *sequence_list){
    //---------------------------------------------------------------------------
    // Feed a string from a file, convert the string into a series of integers
    //---------------------------------------------------------------------------
    clock_t start, end;
    double cpu_time_read, cpu_time_compute, time_malloc;
    start = clock();
    malloc_global_variables();
    end = clock();
    time_malloc += ((double) (end - start)) / CLOCKS_PER_SEC;
    Read *currentRead = malloc(sizeof(Read));
    
    FILE *fp = init_handle_one_file(inputFile);
    // Feed each read and try to detect repeats
    int counter = 0;
    if (numprocs > 1) {
        srand(0);
        for(;;){
            return_one_read(fp, currentRead);
            for (int j = 1; j < numprocs; j++) {
                if (counter == sequence_list[j])
                    srand(j);
            }
            if(currentRead->len == 0) break;
            for(int i=0; i<currentRead->len; i++) {
                orgInputString[i] = currentRead->codedString[i];
            }            
            handle_one_read(currentRead->ID, currentRead->len, tmp_read_cnt, print_alignment);
            counter++;
        }
    } else {
        for(;;){
            start = clock();
            return_one_read(fp, currentRead);
            end = clock();
            cpu_time_read += ((double) (end - start)) / CLOCKS_PER_SEC;
            start = clock();
            if(currentRead->len == 0) break;
            for(int i=0; i<currentRead->len; i++) {
                orgInputString[i] = currentRead->codedString[i];
            }            
            handle_one_read(currentRead->ID, currentRead->len, tmp_read_cnt, print_alignment);
            end = clock();
            cpu_time_compute += ((double) (end - start)) / CLOCKS_PER_SEC;
            counter++;
        }
    }

    fclose(fp);
    fprintf(stderr, "%f\tComputing time\n", cpu_time_compute);
    fprintf(stderr, "%f\tMalloc time\n", time_malloc);
    fprintf(stderr, "%f\tRead time\n", cpu_time_read);
    free_global_variables();
    free(currentRead);
    
    return(tmp_read_cnt);
}
