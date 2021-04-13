#include <vector>
#include <iostream>
#include <fstream>
#include "kdtree.h"
#include <string>
#include <sstream>
#include <chrono>
#include <math.h>
#include <valarray>
#include "Common Tools.h"
// #include "HA_SupportTools.h"


using namespace std;

// Global variables
// How many dimensions
int Dim = 46;
int stepcounter = 0;


// KDTree kd_tree(Dim);

ifstream dataSetFile;
ifstream queryReferenceFile;
std::ofstream timeTable;

std::vector<kdtree::vector_t> recordsCollection;


bool simple_correctness_checker(kdtree::vector_t result, vector<double> test_point, double range, double radius = 0)
{
	result.erase(result.begin());

	if (radius != 0)
	{
		auto ball_x = test_point[0];
		auto ball_y = test_point[1];
		for (auto j = 2; j < test_point.size(); ++j)
		{
			if (j % 2 == 0)
			{
				if (isInsideCircle(ball_x, ball_y, test_point[j], test_point[j + 1], radius))
				{
					auto r = result[j];
					auto t = test_point[j];
					auto comp = abs(r - t);
					if (comp > range)
					{
						return false;
					}
				}
			}
		}
		return true;
	}
	for (int i = 0; i < result.size(); ++i)
	{
		auto r = result[i];
		auto t = test_point[i];
		auto comp = abs(r - t);
		if (comp > range)
		{
			return false;
		}
	}
	return true;
}
void drawVector(vector<double> v, std::string title = " ");
// --- extra
string processing = "";
auto countOfLines = std::count(std::istreambuf_iterator<char>(dataSetFile), std::istreambuf_iterator<char>(), '\n');
std::vector<kdtree::vector_t> referenceForSearch;
// std::ofstream successResults;

// double rangeSet[] = { 5, 6, 7, 8, 9, 10 };

// double rangeSet[] = { 6, 7, 8, 9 };
int qualityThreshold = 10;
double raduisSet[] = { 8, 9, 10, 11, 12, 13, 14 };
double rangeSet[] = { 1, 2, 3, 4, 5 };

std::vector<std::string> recs;

std::vector<kdtree::vector_t> fillTestSet(std::string fileName);
vector<vector<double>> initDataVector_v1(int timeRange, std::string fileName);
//void openFileForSetOfGames(int timeRange = 75)
//{
//	int frame = 0;
//	//int howManyFramesinTree = 0;
//	// =================== Open data set and build tree
//	if (dataSetFile.is_open())
//	{
//		std::cout << "data set file opened" << endl;
//		kdtree::vector_t hPoint;
//		kdtree::vector_t framePoint;
//		hPoint.reserve(Dim);
//		framePoint.reserve(Dim + 1);
//		// Scan data set file and put into kd tree
//		// 105000 = 70min
//		std::string record;
//		while (std::getline(dataSetFile, record))
//		{
//			frame++;
//			if (frame % timeRange == 0)
//			{
//				std::istringstream currentRecord(record);
//				framePoint.push_back(frame);
//				for (auto i = 0; i < Dim; ++i) {
//					// if (i != 0 && i != 1)
//					// {
//					double tempNumber = 0;
//					currentRecord >> tempNumber;
//					hPoint.push_back(tempNumber);
//					framePoint.push_back(tempNumber);
//					// }
//				}
//				try {
//					// insert time
//					// auto startedInsert = std::chrono::high_resolution_clock::now();
//					kd_tree.insert(hPoint, FTagRegionInfo(framePoint));
//					// auto doneInsert = std::chrono::high_resolution_clock::now();
//					// auto elapsedTimeForInsert = std::chrono::duration_cast<std::chrono::microseconds>(doneInsert - startedInsert).count();
//					// std::cout << "Insert time: " << elapsedTimeForInsert << "ms" << endl;
//					stepcounter++;
//					recordsCollection.push_back(hPoint);
//					// howManyFramesinTree++;
//					/*if (stepcounter % 5000 == 0) {
//						processing = "processing tree " + std::to_string(stepcounter) + " of " + std::to_string(countOfLines);
//						std::cout << processing << endl;
//					}*/
//				}
//				catch (kdtree::KeyDuplicateException* e) {} // KeyDuplicateException skip adding to tree
//				hPoint.clear();
//				framePoint.clear();
//			}
//		}
//	}
//	dataSetFile.close();
//}

int main()
{
	// auto dataSetFileName = "G1234.txt"; //"stats_training_set.txt";
	auto dataSetFileName = "d:\\Proccessed Soccer Recordings\\no actions\\Game1-2-3-4\\Game1-2-3-4.txt";
	// auto testSetName = "G1.txt"; //"stats_test_set.txt";
	auto testSetName = "d:\\Proccessed Soccer Recordings\\no actions\\Game5\\game5.txt"; //"stats_test_set.txt";
	timeTable.open("time-table.txt");
	dataSetFile.open(dataSetFileName, ios::in);
	const int frameRate = 25;
	const bool include_ball = true;
	// referenceForSearch = initDataVector_v1(frameRate, testSetName);
	referenceForSearch = initDataVector(frameRate, testSetName, 46, include_ball);

	auto startedBuilding = std::chrono::high_resolution_clock::now(); // TODO: Implement not so clumsy timer

	// openFileForSetOfGames(frameRate);
	auto kd_tree = GetKDtree(frameRate, dataSetFileName, include_ball);
	
	auto doneBuilding = std::chrono::high_resolution_clock::now();
	auto elapsedTimeForBuild = std::chrono::duration_cast<std::chrono::milliseconds>(doneBuilding - startedBuilding).count();

	std::cout << "Build (+Parsing) time: " << elapsedTimeForBuild << "ms" << endl;
	std::cout << "There are " << stepcounter << " elements added into structure" << endl;
	std::vector<pair<float, float>> rangeSuccess;
	rangeSuccess.reserve(sizeof rangeSet / sizeof * rangeSet);

	for (auto radius = 0; radius < sizeof(raduisSet) / sizeof(*raduisSet); ++radius)
	{
		for (auto r2 = 0; r2 < sizeof(rangeSet) / sizeof(*rangeSet); ++r2)
		{
			int rescounter = 0;
			int frameNumber = 0;
			//-------------------- Success statistics gathering
			/*time-Table.open("time_table_Range" +
				std::to_string((int)rangeSet[r2]) +
				"_Radius" +
				std::to_string((int)raduisSet[radius]) +
				".txt");*/
			//---------------------------------------------
			vector<pair<int, int>> frameResultsCounter;

			for (auto i = 0; i < referenceForSearch.size(); ++i) {
				// auto refEx = referenceForSearch[i];
				// auto query = getTestQueryWithGradient(referenceForSearch[i]);
				auto query = getTestQueryAroundBall(referenceForSearch[i], raduisSet[radius], rangeSet[r2]);
				// auto query = getTestQuery(referenceForSearch[i], rangeSet[r2]);

				auto lowerKeys = query.first;
				auto upperKeys = query.second;
				// ========================== Perform query 			
				auto result = kd_tree->range(lowerKeys, upperKeys); // range query
				if (result.size() != 0)
				{
					rescounter++;
					/*int quality = 0;
					for (int nn = 0; nn < lowerKeys.size(); ++nn)
					{
						if (lowerKeys[nn] != minR)
						{
							quality += 1;
						}
					}*/
					//if (quality > qualityThreshold)
					//{
						//timeTable << elapsedTimeForSearch << "," << result.size() << endl;
						// std::cout << result.size() << " results" << endl;
					frameResultsCounter.push_back(make_pair(i, result.size()));
					
					/*  ---- TEMP
					successResults.open("success_findings_" + std::to_string(i)
						+ "_Range_" + radiusString.erase(radiusString.find_first_of('.'), std::string::npos)
						+ "_Radius_" + rangeString.erase(rangeString.find_first_of('.'), std::string::npos)
						+ ".txt");


					successResults << "reference " << endl;
					for (int ref = 0; ref < refEx.size(); ++ref)
					{
						successResults << refEx[ref] << " ";
					}
					successResults << endl;
					successResults << "Query frame: " << i + 1 << " " << std::endl;
					*/
					/*
					successResults << "Range " << std::to_string(fabs(referenceForSearch[i]) - fabs(lowerKeys[i]) << "_Radius" +
						std::to_string((int)raduisSet[radius]) << endl;
					*/

					//for (int k = 0; k < result.size(); ++k)
					//{
					//	if (simple_correctness_checker(result[k].Tag, refEx, rangeSet[r2]/*, raduisSet[radius]*/))
					//	{

					//	}
					//	else {
					//		cout << "There is wrong result" << endl;
					//	}
						/* ---- TEMP
						for (auto j = 0; j < result[0].Tag.size(); ++j)
						{
							// std::cout << result[k].Tag[j] << " ";
							successResults << result[k].Tag[j] << " ";
						}
						successResults << endl;
						*/
						//}
						//} // END of quality condition
				}

				// successResults.close();
				/*
				if (result.size() > 0)
				{
					cout << result.size() << " results" << endl;
					std::cout << "Results:" << endl;
					for (int i = 0; i < result.size(); ++i)
					{
						for (auto j = 0; j < result[0].Tag.size(); ++j)
						{
							std::cout << result[i].Tag[j] << " ";
						}
						std::cout << endl;
					}
				}*/
			}
			if (rescounter > 0)
			{
				// auto radiusString = std::to_string((int)raduisSet[radius]);
				auto rangeString = std::to_string((int)rangeSet[r2]);
				timeTable << "_Range_" + rangeString
					// + "_Radius_" + radiusString
					<< " " << rescounter * 100.0 / referenceForSearch.size() << endl;
				// ========================================================================================================================
				std::cout << "Range = " + rangeString +
					// "; Radius = " + radiusString +
					" " << referenceForSearch.size() << " queries with " << rescounter * 100.0 / referenceForSearch.size() << "% success. "
					<< "With " << rescounter << " results: " << std::endl;
			}

			/*for (auto n = 0; n < frameResultsCounter.size(); ++n)
			{
				std::cout << "frame: " << frameResultsCounter[n].first << " results size " << frameResultsCounter[n].second << endl;
			}*/

			// pair: Range : Success
			// rangeSuccess.push_back(std::make_pair(rangeSet[r2], rescounter * 100.0 / referenceForSearch.size()));
		} // end of Range FOR
	} // end of Radius FOR
	timeTable.close();
	/*for (int i = 0; i < rangeSuccess.size(); ++i)
	{
		cout << "Range: " << rangeSuccess[i].first << " Findings: " << rangeSuccess[i].second << endl;
	}*/
	//dataSetFile.close();
	//dataSetFile.open(dataSetFileName);
	// ---

	// =================== take ref from file	
	//int filecounter = 0;
	//if (queryReferenceFile.is_open())
	//{
	//	stepcounter = 0;
	//	string line;
	//	kdtree::vector_t testPoint;
	//	testPoint.reserve(Dim);
	//	while (getline(queryReferenceFile, line))
	//	{
	//		filecounter++;
	//		stepcounter++;
	//		stringstream ss(line);
	//		for (auto i = 0; i < Dim; ++i)
	//		{
	//			double number;
	//			ss >> number;
	//			testPoint.push_back(number);
	//		}
	//		auto query = getTestQueryAroundBall(testPoint, searchRadius, 2); // getTestQuery(testPoint, 3);
	//		auto lowerKeys = query.first;
	//		auto upperKeys = query.second;

	//		auto startedSearch = std::chrono::high_resolution_clock::now();
	//		auto result = kd_tree.range(lowerKeys, upperKeys); // range query
	//		auto doneSearch = std::chrono::high_resolution_clock::now();


	// manual provided ref value
	/*
	kdtree::vector_t testPoint{ 0.0098, -7.105, 22.645, -0.5243, 7.112, 7.035, 7.854, -3.9984, 7.693, 0.034293, 1.3818, -0.9653, 2.5235, -5.2479, -4.5962, 7.651, 5.3263, -9.219, -1.1466, -2.8812, -0.4641, -9.044, -7.196, 0.024493, -21.196, -0.6566, -6.9678, -11.732, 4.4933, 2.2442, -7.714, -2.1266, 6.2671, 7.672, -8.778, -2.0629, -8.554, 6.1691, 1.1907, -5.2038, -2.5578, -4.8265, -0.0392, -7.028, -0.0735, -4.4541 };
	auto customTest = getTestQueryAroundBall(testPoint, searchRadius, searchRange);
	auto lowerKeys = customTest.first;
	auto upperKeys = customTest.second;
	auto startTestSearch = std::chrono::high_resolution_clock::now();
	auto result = kd_tree.range(lowerKeys, upperKeys); // range query
	auto endTestSearch = std::chrono::high_resolution_clock::now();
	auto elapsedTimeForTestSearch = std::chrono::duration_cast<std::chrono::milliseconds>(endTestSearch - startTestSearch).count();

	cout << "custom test results: " << endl;
	for (int i = 0; i < result.size(); i++)
	{
		for (auto j = 0; j < result[0].Tag.size(); ++j)
		{
			cout << result[i].Tag[j] << " ";
		}
		cout << endl;
	}
	cout << result.size() << " results" << endl;
	*/

	// std::cout << "Test Search time: " << elapsedTimeForTestSearch << "ms" << endl;
	std::cout << "Size of structure " << (sizeof(decltype(recordsCollection.back())) * recordsCollection.size()) / 1024 << "Kb" << endl;

	return 0;
}

std::vector<kdtree::vector_t> fillTestSet(const std::string file_name)
{
	ifstream testSet;
	int frameCounter = 0;
	testSet.open(file_name);
	int howManyTestRecords = 0;
	std::vector<kdtree::vector_t> testSetStructure;
	if (testSet.is_open()) // TODO: Implement universal file opener
	{
		std::string line;
		kdtree::vector_t testPoint;
		testPoint.reserve(Dim);
		int frameRateCounter = 0;
		while (getline(testSet, line))
		{
			frameCounter++;
			stringstream ss(line);
			if (frameRateCounter % 10 == 0)
			{
				for (auto i = 0; i < Dim; ++i)
				{
					double number;
					ss >> number;
					testPoint.push_back(number);
				}
				//howManyTestRecords++;
				testSetStructure.push_back(testPoint);
				testPoint.clear();
			}
		}
	}
	//cout << "!!!!!!    Number of test records is  " << howManyTestRecords << endl;
	testSet.close();
	return testSetStructure;
}

std::vector<kdtree::vector_t> initDataVector_v1(int timeRange = 75, std::string fileName = "")
{
	int frame = 0;
	std::vector<kdtree::vector_t> result_vector;
	ifstream recordFile;
	recordFile.open(fileName, ios::in);
	if (recordFile.is_open())
	{
		std::cout << "opening " << fileName << " file..." << endl;
		std::string record;
		vector<double> framePoint;
		framePoint.reserve(Dim);
		while (std::getline(recordFile, record))
		{
			if (frame % timeRange == 0)
			{
				std::istringstream currentRecord(record);
				for (auto i = 0; i < Dim; ++i) {
					// if (i != 0 && i != 1)
					// {
						double tempNumber = 0;
						currentRecord >> tempNumber;
						framePoint.push_back(tempNumber);
					// }
				}
				result_vector.push_back(framePoint);
				framePoint.clear();
			}
			frame++;
		}
	}
	std::cout << fileName << " extraction complete." << endl;
	recordFile.close();
	return result_vector;
}


