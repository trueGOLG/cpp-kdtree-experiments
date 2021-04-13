#pragma once
#include <vector>
#include <iostream>
#include <sstream>
// #include "Hungarian.h"

using namespace std;
static vector<vector<double>> getDistanceMatrix(vector<double> testRecord, vector<double> dataSetRecord)
{
	auto resultSize = testRecord.size() / 2;
	vector<vector<double>> result(resultSize);
	for (auto& val : result)
	{
		val.resize(resultSize);
	}
	for (int i = 0; i < resultSize; ++i)
	{
		for (int j = 0; j < resultSize; ++j)
		{
			if (i > 10 && j <= 10 ||
				i <= 10 && j > 10)
			{
				result[i][j] = 555;
			}
			else {
				auto x1 = testRecord[2 * i];
				auto y1 = testRecord[2 * i + 1];
				auto x2 = dataSetRecord[2 * j];
				auto y2 = dataSetRecord[2 * j + 1];
				result[i][j] = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
			}
		}
	}
	return result;
}

static vector<double> refine_vector(vector<double>& input_vector, vector<int>& solution)
{
	vector<double> refined(input_vector);
	for (int t = 0; t < solution.size(); ++t)
	{
		refined[t * 2] = input_vector[solution[t] * 2];
		refined[t * 2 + 1] = input_vector[solution[t] * 2 + 1];
	}
	return refined;
}
//void algTask(std::string taskID, vector<vector<double>> querySet, vector<vector<double>> dataSetpart, int range)
//{
//	vector<int> assignment;
//	HungarianAlgorithm HungAlgorithm;
//	int resultCounter3 = 0;
//	int resultCounter4 = 0;
//	int resultCounter5 = 0;
//	int task_size = querySet.size() * dataSetpart.size();
//
//	std::cout << taskID << " started" << endl;
//	int finished = 0;
//
//	for (int i = 0; i < querySet.size(); ++i)
//	{
//		for (auto j = 0; j < dataSetpart.size(); ++j)
//		{
//			auto costMatrix = getDistanceMatrix(querySet[i], dataSetpart[j]);
//			//auto costMatrix = getDistanceMatrix(dataSetpart[j], querySet[i]);
//			double cost = HungAlgorithm.Solve(costMatrix, assignment);
//
//			auto refined = refine_vector(querySet[i], assignment);
//			if (simple_correctness_checker(refined, /*querySet[i]*/ dataSetpart[j], range))
//			{
//				//if (rangeSet[r] == rangeSet[0]) 
//				resultCounter3++;
//				/*if (rangeSet[r] == rangeSet[1]) resultCounter4++;
//				if (rangeSet[r] == rangeSet[2]) resultCounter5++;*/
//			}
//			costMatrix.clear();
//			assignment.clear();
//		}
//	}
//	// cout << "result counter from " << taskID << " is " << resultCounter3 << endl;
//	// cout << "result counter  " << resultCounter3 << "; dataSet part size: " << dataSetpart.size() << endl;
//	const auto dataSetSize = dataSetpart.size();
//	std::cout << taskID << " Finished; For range " << range << " " << resultCounter3 /** 100.0 / dataSetSize*/ << " success items" << endl;
//	/*std::cout << taskID << ": For range " << rangeSet[1] << " " << resultCounter4 * 100.0 / dataSetpart.size() << "% success" << endl;
//	std::cout << taskID << ": For range " << rangeSet[2] << " " << resultCounter5 * 100.0 / dataSetpart.size() << "% success" << endl;*/
//
//}
