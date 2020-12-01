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

char base_int2char(int val){
    switch(val){
        case 0: return('A');
        case 1: return('C');
        case 2: return('G');
        case 3: return('T');
        default: fprintf(stderr, "wrap_around base_int2char: fatal input char %i\n", val); exit(EXIT_FAILURE);
    }
}

int base_char2int(char c){
    switch(c){
        case 'A': return(0);
        case 'C': return(1);
        case 'G': return(2);
        case 'T': return(3);
        default: fprintf(stderr, "wrap_around base_char2int: fatal input char %c\n", c); exit(EXIT_FAILURE);
    }
}

void pretty_print_alignment(char *unit_string, int unit_len, int rep_start, int rep_end, int match_gain, int mismatch_penalty, int indel_penalty){
    
    // tentative
    int MATCH_GAIN          = match_gain;
    int MISMATCH_PENALTY    = mismatch_penalty;
    int INDEL_PENALTY       = indel_penalty;
    
    int rep_unit[MAX_PERIOD];
    int intBase;
    for(int i=0; i<unit_len; i++){   // 1-origin index
        switch(unit_string[i]){
            case 'A': intBase = 0; break;
            case 'C': intBase = 1; break;
            case 'G': intBase = 2; break;
            case 'T': intBase = 3; break;
            default: fprintf(stderr, "base_char2int pretty print: fatal input char %c\n", unit_string[i]);
        }
        rep_unit[i+1] = intBase;    // Shift by 1 for *1*-origin index
    }
    
    int *rep;
    rep = &orgInputString[rep_start-1];
    //rep = &orgInputString[rep_start];
    int rep_len = rep_end - rep_start + 1;
    
    int i, j;
    int next = unit_len+1;
    
    // Initialization
    for(j=0; j<=rep_len; j++){   // Scan rep_unit
        WrapDP[next*0 + j] = 0; // local alignment
    }
    
    int max_wrd = 0;
    int max_i = 0;
    int max_j = 0;
    int val_match, val_mismatch, val_insertion, val_deletion;
    for(i=1; i <= rep_len; i++){        // Scan rep
        for(j=1; j<=unit_len; j++){   // Scan rep_unit
            if( WrapDPsize <= next*i + j ){
                fprintf(stderr, "You need to increse the value of WrapDPsize.\n");
                exit(EXIT_FAILURE);
            }
            if(rep[i] == rep_unit[j]){    // *1*-origin index !!!!
                WrapDP[next*i + j] = WrapDP[next*(i-1) + j-1]  + MATCH_GAIN;
            }else{
                val_mismatch    = WrapDP[next*(i-1) + j-1]  - MISMATCH_PENALTY;
                val_insertion   = WrapDP[next*(i-1) + j]    - INDEL_PENALTY;
                if(j > 1){
                    val_deletion = WrapDP[next*i + j-1] - INDEL_PENALTY;
                    WrapDP[next*i + j] = MAX(0, MAX( MAX( val_mismatch, val_insertion), val_deletion));
                }else{
                    WrapDP[next*i + j] = MAX(0, MAX( val_mismatch, val_insertion));
                }
            }
            if(max_wrd < WrapDP[next*i + j])
            {
                max_wrd = WrapDP[next*i + j];
                max_i = i;
                max_j = j;
            }
        }
        // wrap around
        WrapDP[next*i + 0] = WrapDP[next*i + unit_len];
    }
    
    // trace back the optimal alignment while storing it in the data structure "alignment"
    int Num_matches = 0;
    int Num_mismatches = 0;
    int Num_insertions = 0;
    int Num_deletions  = 0;
    int pos = 0;
    
    i = max_i;
    j = max_j;
    if(j == 0){ j = unit_len; } // 1-origin index
    
    while(i > 0 && WrapDP[next*i + j] > 0){                 // global alignment
        val_match       = WrapDP[next*(i-1) + j-1]  + MATCH_GAIN;
        val_mismatch    = WrapDP[next*(i-1) + j-1]  - MISMATCH_PENALTY;
        val_insertion   = WrapDP[next*(i-1) + j]    - INDEL_PENALTY;
        val_deletion    = WrapDP[next*i + j-1]      - INDEL_PENALTY;
        
        if( max_wrd == val_match && rep[i] == rep_unit[j]){
            alignment_input[pos]   = base_int2char(rep[i]);
            alignment_symbols[pos] = '|';
            alignment_repeats[pos] = base_int2char(rep_unit[j]);
            
            max_wrd -= MATCH_GAIN;
            i--; j--;
            Num_matches++;
            pos++;
        }else if( max_wrd == val_mismatch  && rep[i] != rep_unit[j]){     // mismatch
            alignment_input[pos]   = base_int2char(rep[i]);
            alignment_symbols[pos] = ' ';
            alignment_repeats[pos] = base_int2char(rep_unit[j]);
            
            max_wrd += MISMATCH_PENALTY;
            i--; j--;
            Num_mismatches++;
            pos++;
        }else if( max_wrd == val_deletion){     // deletion
            alignment_input[pos]   = '-';
            alignment_symbols[pos] = ' ';
            alignment_repeats[pos] = base_int2char(rep_unit[j]);

            max_wrd += INDEL_PENALTY;
            j--;
            Num_deletions++;    // Num_insertions++;
            pos++;
        }else if( max_wrd == val_insertion){    // insertion
            alignment_input[pos]   = base_int2char(rep[i]);
            alignment_symbols[pos] = ' ';
            alignment_repeats[pos] = '-';

            max_wrd += INDEL_PENALTY;
            i--;
            Num_insertions++;   // Num_deletions++;
            pos++;
        }else if( max_wrd == 0){
            break;
        }else{
             fprintf(stderr, "fatal error in wrap-around DP max_wrd = %i\n", max_wrd);
             exit(EXIT_FAILURE);
        }
        if(j == 0){
            j = unit_len;
        }
    }

    // print the match gain and mismatch/indel penalties
    printf("match gain = %i, mismatch penalty = %i, indel penalty = %i\n\n", MATCH_GAIN, MISMATCH_PENALTY, INDEL_PENALTY);
    
    // print backward from the end of alignment arrays
    for(int i_start = pos-1; 0 <= i_start; i_start -= ALIGNMENT_WIDTH_PRINTING){
        int i_end;
        if(-1 <= i_start - ALIGNMENT_WIDTH_PRINTING){
            i_end = i_start - ALIGNMENT_WIDTH_PRINTING;
        }else{
            i_end = -1;
        }
        for(int i = i_start; i_end < i; i--){
            printf("%c", alignment_input[i]);
        }
        printf("\n");
        //printf("\t%i-%i\n", rep_start + (pos - 1 - i_start), rep_start + (pos - 1 - i_end - 1));
        
        for(int i = i_start; i_end < i; i--){
            printf("%c", alignment_symbols[i]);
        }
        printf("\n");
        for(int i = i_start; i_end < i; i--){
            printf("%c", alignment_repeats[i]);
        }
        printf("\n\n");
    }
}



// The repeat unit is the input sequence.
// The number of deletions is represented as a positive integer.

// We used a local alignment algorithm in place of a global alignment that outputs lengthy alignments.

void wrap_around_DP_sub( int query_start, int query_end, repeat_in_read *rr, int MATCH_GAIN, int MISMATCH_PENALTY, int INDEL_PENALTY ){
    
    struct timeval s, e;
    gettimeofday(&s, NULL);

    // Initialization
    int unit_len = rr->rep_period;
    
    int rep_unit[MAX_PERIOD];
    for(int i = 0; i < unit_len; i++){
        int intBase;
        switch(rr->string[i]){
            case 'A': intBase = 0; break;
            case 'C': intBase = 1; break;
            case 'G': intBase = 2; break;
            case 'T': intBase = 3; break;
            default: fprintf(stderr, "wrap_around_DP: fatal input char %c at %i (unit len = %i)\n", rr->string[i]), i, rr->rep_period;
        }
        rep_unit[i+1] = intBase;  // Shift by 1 for *1*-origin index
    }
    
    int *rep;
    rep = &orgInputString[query_start];
    int rep_len = query_end - query_start + 1;
    
    int i, j;
    int next = unit_len+1;

    for(j=0; j<=rep_len; j++){   // Scan rep_unit
        WrapDP[next*0 + j] = 0; // local alignment
    }
    
    int max_wrd = 0;
    int max_i = 0;
    int max_j = 0;
    int val_match, val_mismatch, val_insertion, val_deletion;
    for(i=1; i <= rep_len; i++){        // Scan repeat
        for(j=1; j<=unit_len; j++){   // Scan rep_unit
            if( WrapDPsize <= next*i + j ){
                fprintf(stderr, "You need to increse the value of WrapDPsize.\n");
                exit(EXIT_FAILURE);
            }
            if(rep[i] == rep_unit[j]){    // *1*-origin index !!!!
                WrapDP[next*i + j] = WrapDP[next*(i-1) + j-1]  + MATCH_GAIN;
            }else{
                val_mismatch    = WrapDP[next*(i-1) + j-1]  - MISMATCH_PENALTY;
                val_insertion   = WrapDP[next*(i-1) + j]    - INDEL_PENALTY;
                if(j > 1){
                    val_deletion = WrapDP[next*i + j-1] - INDEL_PENALTY;
                    WrapDP[next*i + j] = MAX(0, MAX( MAX( val_mismatch, val_insertion), val_deletion));
                }else{
                    WrapDP[next*i + j] = MAX(0, MAX( val_mismatch, val_insertion));
                }
            }
            if(max_wrd < WrapDP[next*i + j])
            {
                max_wrd = WrapDP[next*i + j];
                max_i = i;
                max_j = j;
            }
        }
        // wrap around
        WrapDP[next*i + 0] = WrapDP[next*i + unit_len];
    }
    
    // trace back the optimal alignment while storing it in the data structure "alignment"
    int Num_matches = 0;
    int Num_mismatches = 0;
    int Num_insertions = 0;
    int Num_deletions  = 0;
    int Num_scanned_unit = 0;
    
    i = max_i;
    j = max_j;
    if(j == 0){ j = unit_len; } // 1-origin index
    
    while(i > 0 && WrapDP[next*i + j] > 0){                 // global alignment
        val_match       = WrapDP[next*(i-1) + j-1]  + MATCH_GAIN;
        val_mismatch    = WrapDP[next*(i-1) + j-1]  - MISMATCH_PENALTY;
        val_insertion   = WrapDP[next*(i-1) + j]    - INDEL_PENALTY;
        val_deletion    = WrapDP[next*i + j-1]      - INDEL_PENALTY;
        
        if( max_wrd == val_match          && rep[i] == rep_unit[j]){
            max_wrd -= MATCH_GAIN;
            i--; j--;
            Num_matches++;
            Num_scanned_unit++;
        }else if( max_wrd == val_mismatch && rep[i] != rep_unit[j]){     // mismatch
            max_wrd += MISMATCH_PENALTY;
            i--; j--;
            Num_mismatches++;
            Num_scanned_unit++;
        }else if( max_wrd == val_deletion){     // deletion
            max_wrd += INDEL_PENALTY;
            j--;
            Num_deletions++;    // Num_insertions++;
            Num_scanned_unit++;
        }else if( max_wrd == val_insertion){    // insertion
            max_wrd += INDEL_PENALTY;
            i--;
            Num_insertions++;
            //Num_scanned_unit++;       // The base of the repeat unit is skipped.
        }else if( max_wrd == 0){
            break;
        }else{
            fprintf(stderr, "fatal error in wrap-around DP max_wrd = %i\n", max_wrd);
            exit(EXIT_FAILURE);
        }
        if(j == 0){
            j = unit_len;
        }
    }
    
    
    // Consider the 1-origin index carefully
    rr->rep_start           = query_start + i + 1;
    rr->rep_end             = query_start + max_i;
    //rr->rep_start           = query_start + i;
    //rr->rep_end             = query_start + max_i - 1;
    rr->repeat_len          = max_i - i; //max_i - 1;
    rr->Num_freq_unit       = (int)Num_scanned_unit/unit_len;
    rr->Num_matches         = Num_matches;
    rr->Num_mismatches      = Num_mismatches;
    rr->Num_insertions      = Num_insertions;
    rr->Num_deletions       = Num_deletions;
    
    rr->match_gain          = MATCH_GAIN;
    rr->mismatch_penalty    = MISMATCH_PENALTY;
    rr->indel_penalty       = INDEL_PENALTY;
 
    gettimeofday(&e, NULL);
    time_wrap_around_DP += (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
}


void wrap_around_DP( int query_start, int query_end, repeat_in_read *rr){
    // Input:   rep_period, string, and string_score
    // Outpur:  repeat_start, repeat_end, Num_freq_unit,
    //          Num_matches, Num_mismatches, Num_insertions, Num_deletions,
    //          match_gain, mismatch_penalty, and indel_penalty
    
    //float rr_ratio = (float)rr->Num_matches / (rr->Num_matches + rr->Num_mismatches + rr->Num_insertions + rr->Num_deletions);
    
    repeat_in_read *tmp_rr, *max_rr;
    tmp_rr = (repeat_in_read*) malloc(sizeof(repeat_in_read));
    if(tmp_rr == NULL){
        fprintf(stderr, "cannot allocate space for tmp_rr.\n");
        exit(EXIT_FAILURE);
    }
    max_rr = (repeat_in_read*) malloc(sizeof(repeat_in_read));
    if(max_rr == NULL){
        fprintf(stderr, "cannot allocate space for max_rr.\n");
        exit(EXIT_FAILURE);
    }
    clear_rr(tmp_rr);
    clear_rr(max_rr);
    float max_ratio = -1;
    float tmp_ratio;
    
    // try MAIN_GAIN = 5, MISMATCH_PENALTY = 1, INDEL_PENALTY = 1
    /*
    set_rr( tmp_rr, rr );
    wrap_around_DP_sub( query_start, query_end, tmp_rr, 5, 1, 1 );
    
    tmp_ratio = (float)tmp_rr->Num_matches / (tmp_rr->Num_matches + tmp_rr->Num_mismatches + tmp_rr->Num_insertions + tmp_rr->Num_deletions);
    if( max_ratio < tmp_ratio ){
        set_rr( max_rr, tmp_rr );
        max_ratio = tmp_ratio;
    }
     */
    
    // try MAIN_GAIN = 1, MISMATCH_PENALTY = 1, INDEL_PENALTY = 3
    set_rr( tmp_rr, rr );
    wrap_around_DP_sub( query_start, query_end, tmp_rr, 1, 1, 3);
    
    tmp_ratio = (float)tmp_rr->Num_matches / (tmp_rr->Num_matches + tmp_rr->Num_mismatches + tmp_rr->Num_insertions + tmp_rr->Num_deletions);
    if(max_ratio < tmp_ratio){
        set_rr( max_rr, tmp_rr );
        max_ratio = tmp_ratio;
    }
    
    // try MAIN_GAIN = 1, MISMATCH_PENALTY = 3, INDEL_PENALTY = 1
    set_rr( tmp_rr, rr );
    wrap_around_DP_sub( query_start, query_end, tmp_rr, 1, 3, 1);
    
    tmp_ratio = (float)tmp_rr->Num_matches / (tmp_rr->Num_matches + tmp_rr->Num_mismatches + tmp_rr->Num_insertions + tmp_rr->Num_deletions);
    if(max_ratio < tmp_ratio){
        set_rr( max_rr, tmp_rr );
        max_ratio = tmp_ratio;
    }
    
    
    set_rr(rr, max_rr);
    
    if(tmp_rr == NULL){
        fprintf(stderr, "cannot free tmp_rr.\n");
        exit(EXIT_FAILURE);
    }else{
        free(tmp_rr);
    }

    if(tmp_rr == NULL){
        fprintf(stderr, "cannot free max_rr.\n");
        exit(EXIT_FAILURE);
    }else{
        free(max_rr);
    }
}
