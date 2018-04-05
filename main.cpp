#include <cstdlib>
#include <cstdio>
#include <limits>
#include <cmath>
#include <omp.h>
#include <time.h>
#include <sys/time.h>
#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <istream>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include <functional>
#include <omp.h>

#define MONTECARLOON false

//#define VERBOSE

#ifdef VVVERBOSE
	#define VVERBOSE
#endif
#ifdef VVERBOSE
	#define VERBOSE
#endif

#define DIR_GRAPH std::string("../graph_checkpointing/")

#define DEBUG(X) std::cout << X << std::endl; std::cout.flush();

#define INF std::numeric_limits<double>::max();
#include "Trace.hpp"

std::vector<double> seeds;

std::vector<std::string> split(const std::string &s, char delim)
{
	std::vector<std::string> tokens;
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss,item,delim))
		tokens.push_back(item);
	return tokens;
}

struct node_t_old
{
	long id;
	double weight;
	std::string label;
        std::vector<std::pair<int,int> > inputs;
        std::vector<std::pair<int,int> > outputs;
        double bottom_level;
	bool checkpoint;
	bool proc;
	bool done;
	bool failed;
};
struct node_t
{
	long id;
	double weight;
	std::string label;
        std::vector<std::pair<int,std::pair<std::string,double> > > inputs;
        std::vector<std::pair<int,std::pair<std::string,double> > > outputs;
	bool checkpoint;
	int proc;
	bool done;
	bool failed;
};

struct graph_t_old
{
	long nb_nodes; //total number of nodes in the graph
        int ntiles;
	std::vector<std::vector<long> > adjacency_list; //two-dimensional array representing the adjacency list
	std::vector<node_t_old> nodes; //the two different possible weights for the tasks
};
struct graph_t
{
	long nb_nodes; //total number of nodes in the graph
	std::vector<node_t> nodes; //the two different possible weights for the tasks
	std::vector<std::vector<long> > schedule; //The list scheduling of the graph
	bool ckpt;
};
std::vector<long> predecessors(graph_t_old*,long);
bool isReady(graph_t_old *G, long task);

bool checkpoint_all(graph_t_old* G, long node)
{
    (void)G; (void)node;
    return true;
}

bool checkpoint_endtasks(graph_t_old* G, long node)
{
    std::string task;
    std::vector<std::string> elems = split(G->nodes[node].label,'_');
    task = elems[0];
    return (task == "POTRF" || task == "TRSM" || task == "GETRF" || task == "TRSMU" || task == "TRSML" || task == "GEQRT" || task == "TSQRT" || task == "UNMQR");
}

double bl_priority(graph_t_old* G, long node)
{
    return G->nodes[node].bottom_level;
}

struct platform_t
{
	int nprocs;//number of the processors
	double lambdaF, lambdaS;//of the fail-stop error and the silent error
	double lambda,f,s;
	double Cd;//cost of checkpoint
	double V;//cost of verification
	double D;//cost of downtime
};

struct pattern_t
{
	double H;//the overhead computed by the formula
	double W;//the time
	int n;//number of segments
	double VG;//the cost of guaranteed verification
	double Cd, Rd;//the cost of the checkpoint,the time of recovery
	double *alpha;//the number of the segment in a w(in this case,there is only one segment)
	double D;//cost of downtime
};

struct simulator_t
{
	long double seedF;
	long double downtime;
	long double start;
	long double horizon;//the stop time

	std::vector<Trace> failstop;
	//std::vector<long> nF;//number of fail-stop error for each processor
	//std::vector<long double *> failstop;//the time that the error strikes
	long global_nF; //number of fail-stop errors for the platform

	double muF;//mu of fail-stop error

};

struct result_t
{
        double time;
        double timeC;
        double nb_faults;
	double H;
	double nRd;
	double nC,nCC; //nb of ckpts
	double k, n;//the number of checkpoint and verification
};

void initSimulator(int nprocs, double lambda, double horizon, simulator_t *s)
{
	s->downtime = 30;
	s->start = 0;
	s->horizon = horizon;

	// Convert lambda to MTBF in year
	if(lambda == 0) s->muF = s->horizon*1000000;
	else s->muF = (1.0/lambda) / ONEYEAR;

	// MTBF for one core
	//s->muF *= p->nprocs;
	//s->muS *= p->nprocs;

	//Initialize the lists of faults for each processor

	s->global_nF = 0;
	for (int i=0; i<nprocs; i++)
	{
		s->failstop.push_back(Trace(s->muF,30,s->horizon));
		s->global_nF += s->failstop[i].nF;
	}
	/*	for (int i=0; i<s->failstop[0].errors.size(); i++)
			std::cout << s->failstop[0].errors[i] << " ";*/

	/*s->nF = std::vector<long>(nprocs,0);
	s->failstop = std::vector<long double *>(nprocs,0);
	for (int i=0; i<nprocs; i++) {
		getFaults(1, rand()+i, s->horizon, s->start, s->muF, s->downtime, &s->failstop[i], &s->nF[i]);
		s->global_nF += s->nF[i];
	}*/
        
}

void nextError(long double cTime, long double *failstop, int *fi, int nF)
{
	// Ensure next fail-stop is after cTime
	while(cTime >= failstop[*fi] && *fi < nF)
		(*fi)++;
}

void readGraph2(graph_t* G, std::string filename, int *n, std::string strategy)
{
	std::map<std::string, long> dict;

//	std::cout << "Reading " << filename << "...\n";
	
	std::ifstream input(filename,std::ios::in);
	std::string line;
	std::getline(input,line);
	while (line.substr(0,7) != "task_id")
		std::getline(input,line);
	std::vector<std::string> split_elems;
	long id=0;
	int nprocs = 0;
	int nC = 0; //number of tasks checkpointed
	while(std::getline(input,line))
	{
		if (line.substr(0,2) == "0,")
			break;
		split_elems = split(line,',');
		dict.insert(std::pair<std::string,long>(split_elems[0],id));
		id++;
		node_t node;
		node.label = split_elems[0];
		node.id = id-1;
		node.weight = std::stod(split_elems[1],nullptr);
		node.proc = std::stoi(split_elems[2],nullptr);
		if (node.proc > nprocs) nprocs = node.proc;
		if (strategy == "ckpt_s1")
			node.checkpoint = (std::stoi(split_elems[3],nullptr)==1);
		else if (strategy == "ckpt_s2")
			node.checkpoint = (std::stoi(split_elems[4],nullptr)==1);
		else if (strategy == "ckpt_all")
			node.checkpoint = (std::stoi(split_elems[5],nullptr)==1);
		else if (strategy == "ckpt_none")
			node.checkpoint = (std::stoi(split_elems[6],nullptr)==1);
		else {
			std::cout << "ERROR: undefined checkpoint strategy. ALL TASKS WILL BE CHECKPOINTED.\n";
			node.checkpoint = true;
		}
		if (node.checkpoint)
			nC++;
		node.done = false; node.failed = false;
		G->nodes.push_back(node);
	}
	G->nb_nodes = id;
	if (strategy == "ckpt_none")
		G->ckpt = false;
	else
		G->ckpt = true;

	std::vector<std::vector<long> > sched = std::vector<std::vector<long> >(G->nb_nodes, std::vector<long>());
	int p=0;
	do
	{
		std::vector<std::string> split_elems = split(line,',');
		for (auto i=1; i<(int)split_elems.size(); i++)
			sched[p].push_back(dict[split_elems[i]]);
		p++;
	} while (std::getline(input,line));
	G->schedule = sched;

	input.clear();
	input.seekg(20,input.beg);
	std::getline(input,line);
	while (line.substr(0,7) != "task_id")
	{
		split_elems = split(line,',');
		for (unsigned i=3; i<split_elems.size(); i++)
		{

			//get file name
			std::size_t firstletter = split_elems[i].find("'") + 1;
			std::size_t lastletter = split_elems[i].find("'",firstletter) - 1;
			std::string edge_name = split_elems[i].substr(firstletter,lastletter-firstletter+1);

			//get file weight
			std::size_t begin = lastletter+4;
			std::size_t end = split_elems[i].find("}",begin);
			double edge_weight = std::stod(split_elems[i].substr(begin,end-begin));

			G->nodes[dict[split_elems[0]]].outputs.push_back(std::pair<long,std::pair<std::string,double> >(dict[split_elems[1]],std::pair<std::string,double>(edge_name,edge_weight) ));
			G->nodes[dict[split_elems[1]]].inputs.push_back(std::pair<long,std::pair<std::string,double> >(dict[split_elems[0]],std::pair<std::string,double>(edge_name,edge_weight) ));
		}
		std::getline(input,line);
	}

	*n = nprocs+1;
	std::cout << nC << " ";

	input.close();

	/*for (auto nn:G->nodes)
	{
		std::cout << nn.id << " " << dict[nn.label] << " " << nn.label << " " << nn.weight << " " << nn.proc << " " << nn.checkpoint << "\n";
		for (auto i:nn.inputs)
			std::cout << i.first << " " << i.second.first << " " << i.second.second << " ";
		std::cout << "\n";
		for (auto i:nn.outputs)
			std::cout << i.first << " " << i.second.first << " " << i.second.second << " ";
		std::cout << "\n";
	}*/

}

std::vector<long> predecessors(graph_t_old* G, long node)
{
	std::vector<long> pred;
	for (unsigned i=0; i<G->nb_nodes; i++)
	{
		if (std::find(G->adjacency_list[i].begin(),G->adjacency_list[i].end(),node) != G->adjacency_list[i].end())
			pred.push_back(i);
	}
	return pred;
}

std::set<long> init_ready(graph_t_old *G)
{
	std::set<long> ready;
	/* OLD
	for (int i=0; i<G->nb_nodes; i++)
	{
		ready.push_back(i);
	}
	for (int i=0; i<G->nb_nodes; i++)
	{
		for (auto neighbor:G->adjacency_list[i])
		{
			for (int j=0; j<ready.size(); j++)
			{
				if (ready[j] == neighbor)
					ready.erase(ready.begin()+j);
			}
		}
	}
	*/
	auto hint = ready.begin();
	for (int i=0; i<G->nb_nodes; i++)
 		if (predecessors(G,i).size() == 0)
			hint = ready.insert(hint,i);
	return ready;
}

bool isReady(graph_t *G, long task)
{
    for (auto p:G->nodes[task].inputs)
    {
        if (!G->nodes[p.first].done)
        {
            return false;
        }
    }
    return true;
}

bool isOnProc(int proc, graph_t* G, long task)
{
	return (proc == G->nodes[task].proc);
}

int startSimulationGraph(int nprocs, simulator_t *s, result_t *r, graph_t *G)
{
	std::vector<int> fi = std::vector<int>(nprocs,0); //index of fail-stop error for each processor
	std::vector<int> si = std::vector<int>(nprocs,0); //index of silent error for each processor

	long exec = 0;//the number of completed tasks
	std::vector<long double> cTime = std::vector<long double>(nprocs,s->start); //The total execution time for each processor

	double VG = 0;

	int nRd = 0;//the number of recoveries
	int nD = 0;//number of downtimes (same as fail-stop errors)
	int nCC = 0; //number of cross-over dependencies (= checkpoints)
	int nC = 0; //number of tasks checkpointed

	//listScheduling will contain the list of the tasks to be executed for each processor
	std::vector<std::vector<long> > listScheduling = G->schedule;
	//currentTask[p]currentTask will contain the index in listScheduling[p] of the task currently being executed (or about to start) for each processor p
	std::vector<int> currentTask = std::vector<int>(nprocs,-1);
	//lastCkpt[p] will contain the index in listScheduling[p] of the last checkpointed task for processor p
	std::vector<int> lastCkpt = std::vector<int>(nprocs,-1);
	//procData[p] will contain the list of already loaded files on proc p
 	std::vector<std::map<std::string,int> > procData = std::vector<std::map<std::string,int> >(nprocs);

	//int test_counter = 0;
	double long global_time = s->start; //application not started
	double long totaltimeC = 0;
	std::vector<long> scheduled = std::vector<long>(nprocs,-1); //List of tasks being executed on each processor
	do {
		#ifdef VERBOSE
		std::cout << "===Time : " << global_time << "===\n";
		#endif
	//	schedule:

		int reset = -1;

                //#pragma omp parallel for
		for (int i=0; i<nprocs; i++) //For all processors, update list of ready tasks if needed
		{
		//DEBUG("reexec" << i);
			//DEBUG(i << " " << cTime[i]);
			if (global_time >= cTime[i]) {
				if (scheduled[i] != -1) { //If the last task scheduled on i is finished
					if (G->nodes[scheduled[i]].failed) //There was an error during the execution
					{
						if (G->ckpt) //in case of ckpt_none strategy, we need to apply a more complicated rollback algorithm as processors are not independent
						{
							//Reset
							exec++; //because listScheduling[i][currentTask[i]] did not ++ yet
							for (auto ind=currentTask[i]; ind>lastCkpt[i]; ind--) {
								G->nodes[listScheduling[i][ind]].done = false;
								G->nodes[listScheduling[i][ind]].failed = false;
								exec--;
							}
							//Rewind to last checkpoint
							currentTask[i] = lastCkpt[i];
						} else {
							#pragma omp critical
							{
							reset = i;
							}
						}
					} else { //Or it was ok, the task finished
						//We mark the task as done, and change the ready tasks
						if (G->nodes[scheduled[i]].done == false)
						{ 
						    G->nodes[scheduled[i]].done = true;
						    exec++;
						    //DEBUG("exec++ " << G->nodes[scheduled[i]].label);
						}
						#ifdef VERBOSE
						std::cout << "task " << G->nodes[scheduled[i]].label << " done.\n";
						#endif
						#ifdef VERBOSE
						std::cout << "Executed : " << exec << "\n";
						#endif
					}
				} else { //no task was scheduled but we may need to reexecute in case of failure during idle time
					if (currentTask[i] > lastCkpt[i]) // if it has already started the execution and CHECK that we have something to redo
					{
						if (G->nodes[listScheduling[i][currentTask[i]]].failed)
						{
							if (G->ckpt) //in case of ckpt_none strategy, we need to apply a more complicated rollback algorithm as processors are not independent
							{
								//Reset
								for (auto ind=currentTask[i]; ind>lastCkpt[i]; ind--) {
									G->nodes[listScheduling[i][ind]].done = false;
									G->nodes[listScheduling[i][ind]].failed = false;
									exec--;
								}
								//Rewind to last checkpoint
								//DEBUG("rollback during idle for proc " << i << " from " << currentTask[i] << " to " << lastCkpt[i]);
								currentTask[i] = lastCkpt[i];
							} else {
								#pragma omp critical
								{
								reset = i;
								}
							}
						}
					}
				}
			}
		}

		if (reset != -1)
		{
			exec = 0;
			for (auto &node: G->nodes) {
				node.failed = false;
				node.done = false;
			}
			for (int i=0; i<nprocs; i++)
			{
				currentTask[i] = -1;
				scheduled[i] = -1;
				lastCkpt[i] = -1;
				procData[i].clear();
				cTime[i] = cTime[reset];
			}
		}

		/*--------------------------------\\
                ||WE SCHEDULE THE TASKS NOW THAT  ||
                ||THE LIST OF READY TASKS HAS BEEN||
                ||UPDATED                         ||
                \\--------------------------------*/


                //#pragma omp parallel for
		for (int i=0; i<nprocs; i++)
		{
		//DEBUG("schedule" << i);
			if (global_time >= cTime[i]) { //Last task scheduled finished, new event
				if (currentTask[i] < (int)listScheduling[i].size()-1) //Last task or not?
				{
					long next_task = listScheduling[i][currentTask[i]+1];
					//Find next task to be scheduled
					if (isReady(G,next_task))
					{
						scheduled[i] = next_task;
						currentTask[i]++;
						cTime[i] = global_time; //set the value of time of proc i to current global time value
					}
					else
						scheduled[i] = -1;
				}
                                else
                                        scheduled[i] = -1;
			}
		}
		#ifdef VVERBOSE
		for (int i=0; i<nprocs; i++)
		{
			if (scheduled[i] != -1)
				std::cout << "P" << i << " has scheduled task " << G->nodes[scheduled[i]].label << "\n";
		}
		#endif
	//	execute:
		double long min_cTime = s->horizon+1;
		// Iterate over processors


                /*------------------------------\\
                ||HERE WE PLAN THE EXECUTION OF ||
                ||THE SCHEDULED TASKS           ||
                \*------------------------------*/

		//#pragma omp parallel for
		for(int i=0; i<nprocs; i++)
		{
		//DEBUG("computation " <<i);
			if(scheduled[i] != -1 && global_time >= cTime[i]) //Check that there is a task to execute and that the task is new
			{
                                double Rtime = 0;
				goto pattern;
				
				disk_recovery:
					#ifdef VERBOSE
					std::cout << "Processor " << i << " entering disk_recovery label at time "<< cTime[i] << "\n";
					#endif
					nRd++;

					///RECOVERY///
					procData[i].clear(); //Remove all loaded files

					Rtime = s->downtime;

					/*WE DON'T DO THE FOLLOWING AS WE WANT TO RECOVER AS LATE AS POSSIBLE (==> READ INPUT part) THE FILES

					//Computation of recovery time: we recover all data used by tasks after last checkpointed from SAME processor
					//because if data from another processor is needed it will be read from disk at the execution (see below)
					for (auto ind=lastCkpt[i]+1; ind<(int)listScheduling[i].size(); ind++)
					{
						for (auto data:G->nodes[listScheduling[i][ind]].inputs)
						{
							for (auto ind2=0; ind2<=lastCkpt[i]; ind2++)
							{
								if (data.first == listScheduling[i][ind2]) //the tasks correspond (i.e. we have an edge crossing the checkpoint)
									Rtime += data.second;
							}
						}
					}*/

					Rtime /= ONEYEAR;

					//Rtime is only downtime so no error during recovery
/*
					nextError(cTime[i], s->failstop[i], &fi[i], s->nF[i]);

					if(s->failstop[i][fi[i]] - cTime[i] < Rtime)
					{
                                            #ifdef VVERBOSE
						std::cout << "Fail-stop error during recover after execution of " << G->nodes[scheduled[i]].label << " at time " << s->failstop[i][fi[i]] << "\n";
					    #endif
                                            cTime[i] += (s->failstop[i][fi[i]] - cTime[i]);
						goto disk_recovery;
					}
					else*/
						cTime[i] += Rtime;
					#ifdef VVERBOSE
					std::cout << "Recovered after " << G->nodes[scheduled[i]].label << " => time " << cTime[i] << "\n";
					#endif
					G->nodes[scheduled[i]].failed = true; //Flag to know that we need to reschedule some tasks at next step
					continue; //Go to next processor
				pattern:
					#ifdef VVERBOSE
					//std::cout << "Processor " << i << " entering pattern label\n";
					#endif
					//nextError(cTime[i], s->failstop[i], &fi[i], s->nF[i]);
					double nextError = s->failstop[i].next(cTime[i]);

					double VV = VG;
					double weight = G->nodes[scheduled[i]].weight;

					///READ INPUT///
					//we add the cost of all task inputs until next checkpoint to simulate monte-carlo method
					if (MONTECARLOON && currentTask[i]-1 == lastCkpt[i]) //Just after a checkpoint/recovery
					{
						std::map<std::string,int> readData;
						for (auto data:G->nodes[scheduled[i]].inputs)
						{
							auto ret = procData[i].insert(std::make_pair(data.second.first,1));
							//We check all input files and insert them in the list in they are not present (+ pay the cost)
							if (ret.second) {
								readData.insert(std::make_pair(data.second.first,1));
								weight += data.second.second;
								#ifdef VERBOSE
									std::cout << "Data " << data.second.first << " was read (cost " << data.second.second << ").\n";
								#endif
							}
						}
						for (auto data:G->nodes[scheduled[i]].outputs)
							readData.insert(std::make_pair(data.second.first,1));
						for (auto ind=currentTask[i]; true; ind++)
						{
							if (G->nodes[listScheduling[i][ind]].checkpoint) //stop when ckpt found
								break;
							if (listScheduling[i].size() <= (unsigned)ind+1) //security check if reached the end
								break;
							for (auto data:G->nodes[listScheduling[i][ind+1]].inputs)
							{
								auto ret = readData.insert(std::make_pair(data.second.first,1)); //Check in READDATA not PROCDATA
								if (ret.second)
								{
									procData[i].insert(std::make_pair(data.second.first,1));
									weight += data.second.second;
									#ifdef VERBOSE
										std::cout << "Data " << data.second.first << " was read (cost " << data.second.second << ").\n";
									#endif
								}
							}
							for (auto data:G->nodes[listScheduling[i][ind+1]].outputs)
								readData.insert(std::make_pair(data.second.first,1));
						}
					} else { //CLASSICAL READ AS LATE AS POSSIBLE
						for (auto data:G->nodes[scheduled[i]].inputs)
						{
							auto ret = procData[i].insert(std::make_pair(data.second.first,1));
							//We check all input files and insert them in the list in they are not present (+ pay the cost)
							if (ret.second) {
								weight += data.second.second;
								#ifdef VERBOSE
									std::cout << "Data " << data.second.first << " was read (cost " << data.second.second << ").\n";
								#endif
							}
						}
					}
					weight /= ONEYEAR;

					// If there is a fail-stop error
					if(nextError - cTime[i] < weight + VV)
					{
						#ifdef VVERBOSE
						std::cout << "Fail-stop error during execution of " << G->nodes[scheduled[i]].label << " at time " << nextError << "\n";
						#endif
						nD++;
						cTime[i] += (nextError - cTime[i]);
						goto disk_recovery;
					}
					// Else if there is a silent error
                                        // NO SILENT ERROR FOR NOW
					/*else if(s->silent[i][si[i]] - cTime[i] < weight)//
					{
						#ifdef VVERBOSE
						std::cout << "Silent error during execution of " << G->nodes[scheduled[i]].label << " at time " << s->silent[i][si[i]] << "\n";
						#endif
						cTime[i] += (weight + VV);
						sn++;
						goto disk_recovery;
					}*/
					// Else if no error
					else
					{
						cTime[i] += (weight + VV);

						///DATA
						for (auto data:G->nodes[scheduled[i]].outputs)
							procData[i].insert(std::make_pair(data.second.first,1)); //We have the outputs file in memory
                                                //DEBUG("exec++ task " << G->nodes[scheduled[i]].label << " processor " << i);
						//exec++;

					}

					//DEBUG("ckpt " <<i);

                                        //Checkpoint des données de sortie
					///CHECKPOINT///
					double ckptTime = 0;
					std::map<std::string,int> ckptData;
					//We first checkpoint the cross-over dependencies
					if (G->ckpt) { //FOR CKPT NONE STRATEGY
					for (auto data:G->nodes[listScheduling[i][currentTask[i]]].outputs)
					{
						if (! isOnProc(i,G,data.first) || G->nodes[scheduled[i]].checkpoint) //Checkpoint cross-over (and others if checkpoint taken)
						{
							//check if the file was already checkpointed
							auto ret = ckptData.insert( std::make_pair(data.second.first,1) );
							if (ret.second)
							{
								ckptTime += data.second.second;
								#ifdef VERBOSE
									std::cout << "Data " << data.second.first << " was checkpointed (cost " << data.second.second << ").\n";
								#endif
							}
						}
					}

					if (ckptTime > 0) nCC++;

				    	if(G->nodes[scheduled[i]].checkpoint) //if the task has to be checkpointed, we add the dependencies on same processor
				    	{
						//Computation of checkpointing cost
						//We start from last checkpoint up to task before current task to checkpoint their data
						for (auto ind=lastCkpt[i]+1; ind < currentTask[i]; ind++)
						{
							for (auto data:G->nodes[listScheduling[i][ind]].outputs)
							{
								//we look if data will be needed in next tasks
								for (auto ind2=currentTask[i]+1; ind2 < (int)listScheduling[i].size(); ind2++)
								{
									//we check if a task that will be executed after (listScheduling[i][ind2]) will need the data from a not checkpointed task (listScheduling[i][ind])
									if (listScheduling[i][ind2] == data.first)
									{
										//check if the file was already checkpointed
										auto ret = ckptData.insert( std::make_pair(data.second.first,1) );
										if (ret.second)
										{
											ckptTime += data.second.second;
											#ifdef VERBOSE
												std::cout << "Data " << data.second.first << " was checkpointed (cost " << data.second.second << ").\n";
											#endif
										}
									}
								}
							}
						}
						//And we always checkpoint every output data for the task being checkpointed => done 28 lines above
						nC++;
                                        }
					}
					ckptTime /= ONEYEAR;

					// Perform disk checkpoint
					if(nextError - cTime[i] < ckptTime)//if an error occur
					{
						#ifdef VVERBOSE
						std::cout << "Fail-stop error during checkpoint of " << G->nodes[scheduled[i]].label << " at time " << nextError << "\n";
						#endif
						cTime[i] += (nextError - cTime[i]);
						totaltimeC += (nextError - cTime[i]);
						/*exec--; //The task need to be redone because of unexpected failure...
                                                        DEBUG("exec-- " << G->nodes[scheduled[i]].label);*/
                                                        goto disk_recovery;
					}
					else
					{
						// We reached the final disk checkpoint of the pattern without error!
						cTime[i] += ckptTime; ///Change value according to input data
						totaltimeC += ckptTime;
						if (G->nodes[scheduled[i]].checkpoint)
						{
							lastCkpt[i] = currentTask[i]; //Update the last checkpointed task index
							//We clear the data loaded after every checkpoint
							procData[i].clear();
						}

					}
			} else {
				//CHECK ERRORS BETWEEN cTime[i] and globaltime
				/*double nextError = s->failstop[i].next(cTime[i]);
				while (nextError <= global_time) //in case several errors strike during idle
				{
					//DEBUG("failstop during idle for proc " << i << " after task " << currentTask[i] << " (last ckpt: " << lastCkpt[i] << ")");
					nRd++;
					procData[i].clear();
					cTime[i] = nextError + s->downtime/ONEYEAR;
					if (currentTask[i] > lastCkpt[i]) //if some work was lost
						G->nodes[listScheduling[i][currentTask[i]]].failed = true; //Flag to know that we need to reschedule some tasks at next step
					nextError = s->failstop[i].next(cTime[i]);
				}*/
			}
		}

		#ifdef VVVERBOSE
		for (int i=0; i<nprocs; i++)
		{
			std::cout << "Proc " << i << ": \t";
			for (auto t:listScheduling[i])
			{
				if (t == listScheduling[i][currentTask[i]])
					std::cout << "!";
				std::cout << G->nodes[t].label << " ";
			}
			std::cout << "\n";
		}
		#endif


                if (exec < G->nb_nodes)
                {
		for (int i=0; i<nprocs; i++) //We get the next value of time for which some task could be scheduled
		{
			#ifdef VVVERBOSE
			  DEBUG("processor " << i << " : time " << cTime[i] << ", scheduled: " << scheduled[i]);
			#endif
			if (cTime[i] < min_cTime && scheduled[i] != -1)
				min_cTime = cTime[i];
		}
		global_time = min_cTime;
                }
		//if (global_time >= 4)
		//	std::cerr << "STOOOOOOOP";
		//test_counter++;
		//DEBUG(exec << " <? " << G->nb_nodes);

	}
	while(exec < G->nb_nodes && global_time < s->horizon);
        #ifdef VERBOSE
           std::cout << "End time : " << global_time << "\n";
        #endif
	if(G->nb_nodes - exec > 0)
	{
		//printf("ERROR: application took too long. Increase Horizon.\n");
		r->H = -1;
		r->nb_faults = nRd;
	}
	else
	{
		// Done!
		double totaltime = (global_time - s->start);
                totaltime *= ONEYEAR;
                totaltime /= ONEHOUR;
		double totalworktime = 0;
		for (int i=0; i<G->nb_nodes; i++)
		{
			totalworktime += G->nodes[i].weight;
		}
		totalworktime /= ONEHOUR;
                r->time = totaltime;
		r->timeC = totaltimeC * ONEYEAR / ONEHOUR;
                r->nb_faults = nRd;
                //DEBUG(totaltime << " " << totalworktime);
		// Normalized overhead with respect to the total worktime
		r->H = totaltime / totalworktime;
		// Number of errors per year
		r->nRd = (double)nRd / totaltime;
		r->nC = nC;
		r->nCC = nCC;
	}
	return 0;
}

void init()
{
	struct timeval time;
	gettimeofday(&time, NULL);
	srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
}

void printParam(std::string filename,double* la)
{
	std::vector<std::string> split_elems = split(filename,'_');
	std::cout << split_elems[0].substr(5) << " " << split_elems[1] << " " << split_elems[2].substr(4) << " " << split_elems[3].substr(3) << " " << split_elems[4].substr(3) << " ";

	*la = std::stod(split_elems[3].substr(3));

}

void printResults(result_t *r)
{
	/*std::cout << "------RESULTS------\n";
	std::cout << "Total time: " << r->time << "\n";
	std::cout << "Number of fail-stop errors: " << r->nb_faults << "\n";
	std::cout << "Overhead: " << r->H << "\n";*/
	if (r->H == -1)
	{
		r->time = -1;
	}
	std::cout << r->time << " " << r->nb_faults << " " << r->timeC << " " << r->nC << " " << r->nCC << " ";
}

int main(int argc, char **argv)
{

	if (argc != 4)
	{
		std::cout << "Wrong number of parameters. Usage: " << argv[0] << " input_file [ckpt_s1|ckpt_s2|ckpt_all|ckpt_none|every|almost] horizon.\n";
		return 1;
	}

	init();

	int nprocs = 0;
	double lambda = 0;
	double horizon = atof(argv[3]);
	result_t result;

	std::string test_file(argv[1]);
	printParam(test_file,&lambda);
	//lambda = 0.00001;

	std::string mode(argv[2]);
	std::string strat;

	if (mode != "every" && mode != "almost")
		strat = mode;
	else
		strat = "ckpt_s1";

	//New variables
	beginSimulation:
	graph_t test_graph;
	simulator_t simulator;

	//Read, run, print, clean
	readGraph2(&test_graph,test_file,&nprocs,strat);
	initSimulator(nprocs,lambda,horizon,&simulator);
	startSimulationGraph(nprocs,&simulator,&result,&test_graph);
	printResults(&result);

	if ((mode == "every" || mode == "almost") && strat == "ckpt_s1")
	{
		strat = "ckpt_s2";
		goto beginSimulation;
	} else if ((mode == "every" || mode == "almost") && strat == "ckpt_s2") {
		strat = "ckpt_all";
		goto beginSimulation;
	} else if (mode == "every" && strat == "ckpt_all") {
		strat = "ckpt_none";
		goto beginSimulation;
	}

	std::cout << std::endl;

	return 0;
}
