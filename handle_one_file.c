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

void read_sequences(FILE *fp, int *seq_lens, long int *seq_start_pointers, int *total_seqs_len, int myid, int numprocs){
    // Return 0 if it is the last read
    char *s = (char *)malloc(sizeof(char)*BLK);
    int len;
    int cnt=0;
    int string_cnt=0;
    int total_string_cnt=0;
    long int seq_pointer = 0;

    while (fgets(s, BLK, fp) != NULL) { // Feed a string of size BLK from fp into string s
        if(s[0] == '>'){
            seq_start_pointers[cnt] = seq_pointer;
            seq_lens[cnt] = string_cnt;
            total_string_cnt += string_cnt;
            cnt++;
            string_cnt=0;
        } else {
            // Feed the string
            for(int i=0; s[i]!='\0' && s[i]!='\n' && s[i]!='\r'; i++){
            //for(i=0; s[i]!='\0' && s[i]!='\n'; i++){
                string_cnt++;
                if( MAX_INPUT_LENGTH <= string_cnt){
                    fprintf(stderr, "fatal error: The length %d is tentatively at most %i.\nread ID = %d\nSet MAX_INPUT_LENGTH to a larger value", string_cnt, MAX_INPUT_LENGTH, cnt);
                    free_global_variables_and_exit();
                    exit(EXIT_FAILURE);
                }
            }
        }
        seq_pointer = ftell(fp);
    }

    seq_lens[cnt] = string_cnt;
    *total_seqs_len = total_string_cnt;
    rewind(fp);
}

int return_all_reads(FILE *fp, Read *readList, int myid, int numprocs){
    // Return 0 if it is the last read
    char *s = (char *)malloc(sizeof(char)*BLK);
    int i;
    int cnt=-1;
    int string_cnt=0;
    
    while (fgets(s, BLK, fp) != NULL) { // Feed a string of size BLK from fp into string s
        if(s[0] == '>'){
            cnt++;
            string_cnt=0;
            for(i=1; s[i]!='\0' && s[i]!='\n' && s[i]!='\r' && i<BLK; i++) {
                readList[cnt].ID[i-1] = s[i];
            }
            
            readList[cnt].ID[i-1] = '\0';
        } else {
            // Feed the string
            for(i=0; s[i]!='\0' && s[i]!='\n' && s[i]!='\r'; i++){
            //for(i=0; s[i]!='\0' && s[i]!='\n'; i++){
                readList[cnt].codedString[string_cnt++] = char2int(s[i]);
                if( MAX_INPUT_LENGTH <= string_cnt){
                    fprintf(stderr, "fatal error: The length %d is tentatively at most %i.\nread ID = %s\nSet MAX_INPUT_LENGTH to a larger value", string_cnt, MAX_INPUT_LENGTH, readList[cnt].ID);
                    free_global_variables_and_exit();
                    exit(EXIT_FAILURE);
                }
            }
            readList[cnt].len = string_cnt;
        }
        
    }

    return cnt+1;
}


int handle_one_file(char *inputFile, int print_alignment, int myid, int numprocs, int print_computation_time){
    //---------------------------------------------------------------------------
    // Feed a string from a file, convert the string into a series of integers
    //---------------------------------------------------------------------------
    double handle_read_time, communication_time, program_time, es_time, malloc_time, loadbalance_time;
    double start, end, all_start, all_end;

    if (print_computation_time) {
        all_start = MPI_Wtime();
        start = MPI_Wtime();
    }

    malloc_global_variables();
    int num_sequences;
    int *start_sequences, *end_sequences, *seq_lens;
    long int *out_ptrs = (long int *) malloc(numprocs * sizeof(long int));
    long int *seq_start_pointers;
    char *s = (char *)malloc(sizeof(char)*BLK);

    if (print_computation_time) {
        end = MPI_Wtime();
        malloc_time += end-start;
        start = MPI_Wtime();
    }
    FILE *fp = init_handle_one_file(inputFile);
    if(myid == 0){
        // Get the number of sequences in sequences file
        while (fgets(s, sizeof(s), fp)){
            if (s[0] == '>'){
                num_sequences++;
            }
        }
        fprintf(stderr, "number of sequences %d\n", num_sequences);
        rewind(fp);
    }
    if (print_computation_time) {
        end = MPI_Wtime();
        es_time += end-start;
        start = MPI_Wtime();
    }
    MPI_Bcast(
        &num_sequences,
        1,
        MPI_INT,
        0,
        MPI_COMM_WORLD);

    if (print_computation_time) {
        end = MPI_Wtime();
        communication_time += end-start;
        start = MPI_Wtime();
    }

    // Sequences lengths
    seq_lens = (int *) malloc(num_sequences * sizeof(int));
    // Sequences start positions in sequences file
    seq_start_pointers = (long int *) malloc(num_sequences * sizeof(long int));

    // Start sequences for each process
    start_sequences = (int *) malloc(numprocs * sizeof(int));
    // End sequences for each process
    end_sequences = (int *) malloc(numprocs * sizeof(int));

    if (print_computation_time) {
        end = MPI_Wtime();
        malloc_time += end-start;
        start = MPI_Wtime();
    }

    // Feed each read and try to detect repeats
    if (myid == 0) {
        int total_seqs_len;
        read_sequences(fp, seq_lens, seq_start_pointers, &total_seqs_len, myid, numprocs);
        long int seq_lens_per_proc = total_seqs_len / numprocs;
        int current_proc = 0;
        long int len_acu = 0;
        
        for (int seq = 0; seq < num_sequences; seq++) {
            if (len_acu >= seq_lens_per_proc * current_proc) {
                if (current_proc == 0) {
                    end_sequences[numprocs-1] = num_sequences;
                } else {
                    end_sequences[current_proc-1] = seq;
                }
                start_sequences[current_proc] = seq;
                current_proc++;
            }
            len_acu += seq_lens[seq];
        }
    }

    if (print_computation_time) {
        end = MPI_Wtime();
        loadbalance_time += end-start;
        start = MPI_Wtime();
    }
    // Broadcast sequences lengths from process 0 to all processes
    MPI_Bcast(seq_lens, num_sequences, MPI_INT, 0, MPI_COMM_WORLD);
    // Broadcast sequences start positions in sequences file from process 0 to all processes
    MPI_Bcast(seq_start_pointers, num_sequences, MPI_LONG, 0, MPI_COMM_WORLD);
    // Broadcast start sequences from process 0 to all processes
    MPI_Bcast(start_sequences, numprocs, MPI_INT, 0, MPI_COMM_WORLD);
    // Broadcast end sequences from process 0 to all processes
    MPI_Bcast(end_sequences, numprocs, MPI_INT, 0, MPI_COMM_WORLD);
    if (print_computation_time) {
        end = MPI_Wtime();
        communication_time += end-start;
    }

    // Get sequences to do for every process
    int num_start_seq = start_sequences[myid];
    int num_end_seq = end_sequences[myid];
    int current_seq = num_start_seq;

    
    // Move process to each initial position in file
    fseek(fp, seq_start_pointers[num_start_seq], SEEK_SET);
    memset(s, '\0', sizeof(s));
    int i;
    int cnt=0;
    int string_cnt=0;

    char output_filename[4096] = "";
    sprintf(output_filename, "result-%d.txt", myid);
    FILE *output_file = fopen(output_filename, "w+");

    srand(myid);

    Read currentRead;
    if (print_computation_time) {
        start = MPI_Wtime();
    }
    //fprintf(stderr, "initial seq for process %d: %d\n", myid, current_seq);
    while (fgets(s, BLK, fp) != NULL && (current_seq < num_end_seq)) { // Feed a string of size BLK from fp into string s
        if(s[0] == '>'){
            if (cnt > 0 ) {
                for(int j=0; j < currentRead.len; j++) {
                    orgInputString[j] = currentRead.codedString[j];
                }
                handle_one_read(currentRead.ID, currentRead.len, num_sequences, print_alignment, myid, output_file);
                current_seq++;
            }
            cnt++;
            string_cnt=0;
            for(i=1; s[i]!='\0' && s[i]!='\n' && s[i]!='\r' && i<BLK; i++) {
                currentRead.ID[i-1] = s[i];
            }
            currentRead.ID[i-1] = '\0';
        } else {
            // Feed the string
            for(i=0; s[i]!='\0' && s[i]!='\n' && s[i]!='\r'; i++){
                currentRead.codedString[string_cnt++] = char2int(s[i]);
                if( MAX_INPUT_LENGTH <= string_cnt){
                    fprintf(stderr, "fatal error: The length %d is tentatively at most %i.\nread ID = %s\nSet MAX_INPUT_LENGTH to a larger value", string_cnt, MAX_INPUT_LENGTH, currentRead.ID);
                    free_global_variables_and_exit();
                    exit(EXIT_FAILURE);
                }
            }
            currentRead.len = string_cnt;

        }
    }

    if (myid == numprocs - 1) {
        for(int j=0; j < currentRead.len; j++) {
            orgInputString[j] = currentRead.codedString[j];
        }
        handle_one_read(currentRead.ID, currentRead.len, num_sequences, print_alignment, myid, output_file);
        current_seq++;
    }

    if (print_computation_time) {
        end = MPI_Wtime();
        handle_read_time = end-start;
    }

    fclose(fp);
    long int out_seq_ptr = 0;
    int current_proc = 0;
    char *result_filename = "result.txt";
    out_ptrs = (long int *) malloc(numprocs * sizeof(long int));
    MPI_Barrier(MPI_COMM_WORLD);

    long int seq_pointer = ftell(output_file);
    rewind(output_file);

    if (print_computation_time) {
        start = MPI_Wtime();
    }

    MPI_Allgather(&seq_pointer, 1, MPI_LONG, out_ptrs, 1, MPI_LONG, MPI_COMM_WORLD);

    if (print_computation_time) {
        end = MPI_Wtime();
        communication_time += end-start;
    }
    for(current_proc = 0; current_proc < myid; current_proc++){
        out_seq_ptr += out_ptrs[current_proc];
    }
    
    if (print_computation_time) {
        start = MPI_Wtime();
    }
    MPI_File result_file;
    MPI_File_open(MPI_COMM_WORLD, result_filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &result_file);
    MPI_File_seek(result_file, out_seq_ptr, MPI_SEEK_SET);
    while(fgets(s, BLK, output_file) != NULL) {
        MPI_File_write(result_file, s, strlen(s), MPI_CHAR, MPI_STATUS_IGNORE);
    }
    MPI_File_close(&result_file);
    if (print_computation_time) {
        end = MPI_Wtime();
        es_time += end-start;
        all_end = MPI_Wtime();
        program_time = all_end - all_start;
    }

    if (print_computation_time) {
        double avgtime_comm, avgtime_read, avgtime_prg, avgtime_es, avgtime_malloc, avgtime_lb;
        start = MPI_Wtime();
        MPI_Reduce(&communication_time, &avgtime_comm, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
        MPI_Reduce(&es_time, &avgtime_es, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
        MPI_Reduce(&malloc_time, &avgtime_malloc, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
        MPI_Reduce(&loadbalance_time, &avgtime_lb, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
        MPI_Reduce(&handle_read_time, &avgtime_read, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
        MPI_Reduce(&program_time, &avgtime_prg, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
        end = MPI_Wtime();
        
        if (myid==0) {
            avgtime_comm += (end-start)/numprocs;
            fprintf(stderr, "avg. time communications: %f\n", avgtime_comm/numprocs);
            fprintf(stderr, "avg. time handle read: %f\n", avgtime_read/numprocs);
            fprintf(stderr, "avg. time malloc: %f\n", avgtime_malloc/numprocs);
            fprintf(stderr, "avg. time loadbalance: %f\n", avgtime_lb/numprocs);
            fprintf(stderr, "avg. time handle I/O: %f\n", avgtime_es/numprocs);
            fprintf(stderr, "avg. time handle_one_file function: %f\n", avgtime_prg/numprocs);
        }

    }

    fclose(output_file);
    free_global_variables();
    free(out_ptrs);
    free(seq_start_pointers);
    free(start_sequences);
    free(end_sequences);
    free(seq_lens);
    return(tmp_read_cnt);
}