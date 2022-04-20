//
//  count_match.cpp
//
//
//  Created by Kazuki Ichikawa
//
#include <random>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <strstream>
#include <vector>
#include <sstream>

using namespace std;

mt19937 mt;

string rand_string()
{
	string tmp_s;

	int tmp = (int)(mt()%4);
	
	if(tmp == 0)
	{
		tmp_s = "A";
	}
	else if(tmp == 1)
	{
		tmp_s = "T";
	}
	else if(tmp == 2)
	{
		tmp_s = "C";
	}
	else
	{
		tmp_s = "G";
	}

	
	return tmp_s;
}

void rand_fasta(const char* out_name, const char* rep_name, int rep_length, int block, double mis_rate, double ins_rate, double del_rate, int pre, int post, int loop)
{
	random_device rnd;
    mt.seed(rnd()); 
	
	int rep_len = rep_length * block;
	
	int mis_n = (int)round(((double)(rep_length * block * mis_rate)/100));
	int ins_n = (int)round(((double)(rep_length * block * ins_rate)/100));
	int del_n = (int)round(((double)(rep_length * block * del_rate)/100));
	
	ofstream ofs;
	ofs.open(out_name);
	ofstream ofs2;
	ofs2.open(rep_name);
	int random;

	for(int L=0; L<loop; L++)
	{
		ofs << ">"  << L << endl;
		
		random  = (rand() % (post + 1 - pre)) + pre;
		for(int i=0; i<random; i++)
		{
			ofs << rand_string();
		}
		
		int N = random + 1;
		
		vector<int> error_row;
		
		for(int i=0; i<rep_len; i++)
		{
			error_row.push_back(0);
		}
		
		for(int i=0; i<mis_n; i++)
		{
			while(1)
			{
				int tmp_n = (int)(mt() % rep_len);
				
				if(error_row[tmp_n] == 0)
				{
					error_row[tmp_n] = 1;
					break;
				}
			}
		}

		for(int i=0; i<ins_n; i++)
		{
			while(1)
			{
				int tmp_n = (int)(mt() % rep_len);

				if(error_row[tmp_n] == 0)
				{
					error_row[tmp_n] = 2;
					break;
				}
			}
		}
		
		for(int i=0; i<del_n; i++)
		{
			while(1)
			{
				int tmp_n = (int)(mt() % rep_len);

				if(error_row[tmp_n] == 0)
				{
					error_row[tmp_n] = 3;
					break;
				}
			}
		}

		string line;
		
		while(1)
		{
			line = "";
			
			for(int i=0; i<rep_length; i++)
			{
				line += rand_string();
			}

			bool all_dif = true;
			
			for(int i=1; i<rep_length; i++)
			{
				if(rep_length % i == 0)
				{
					bool tmp_same = true;
					
					int div = rep_length / i;
					
					string sub_s = line.substr(0, i);
					
					for(int j=1; j<div; j++)
					{
						string tmp_s = line.substr(j * i, i);
						
						if(sub_s != tmp_s)
						{
							tmp_same = false;
							break;
						}
					}
					
					if(tmp_same)
					{
						all_dif = false;
						break;
					}
				}
			}
			
			if(all_dif)
			{
				break;
			}
		}
				
		ofs2 << line << endl;

		int rep_total = 0;
		
		for(int i=0; i<block; i++)
		{
			string tmp = line;
			
			for(int j=0; j<(int)tmp.size(); j++)
			{
				if(error_row[rep_total] == 1)
				{
					string mis;
					while(1)
					{
						mis = rand_string();
						string test;
						test = line[j];
						if(mis != test)
						break;
					}
					ofs << mis;
					N++;
				}
				else if(error_row[rep_total] == 2)
				{
					ofs << tmp[j];
					ofs << rand_string();
					N++;
					N++;
				}
				else if(error_row[rep_total] == 3)
				{
				}
				else
				{
					ofs << tmp[j];
					N++;
				}
				rep_total++;
			}
		}
		
		for(int i=0; i<random; i++)
		{
			ofs << rand_string();
		}
		
		ofs << endl;
	}
}

int main(int argc, char** argv)
{
	rand_fasta(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]), atoi(argv[8]), atoi(argv[9]), atoi(argv[10]));
}