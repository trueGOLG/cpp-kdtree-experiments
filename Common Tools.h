#pragma once
#include "Hungarian.h"
#include "kdtree.h"
#include "Common Tools.h"
#include <vector>
#include <algorithm>

static std::vector<kdtree::vector_t> initDataVector(int timeRange = 75, std::string fileName = "", int Dim = 46, bool include_Ball = true)
{
	int frame = 0;
	std::vector<kdtree::vector_t> result_vector;
	ifstream recordFile;
	recordFile.open(fileName);
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
					double tempNumber = 0;
					currentRecord >> tempNumber;
					framePoint.push_back(tempNumber);
				}
				// remove ball coordinates
				if (!include_Ball)
				{
					std::rotate(framePoint.begin(), framePoint.begin() + 2, framePoint.end());
					framePoint.pop_back(); framePoint.pop_back();
				}
				// -----------------------
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
static std::vector<kdtree::vector_t> initDataVectorWithCentroids(vector<double> centroids, int timeRange = 75, std::string fileName = "", int Dim = 46, bool include_Ball = true)
{
	HungarianAlgorithm HungAlg;
	int frame = 0;
	std::vector<kdtree::vector_t> result_vector;
	ifstream recordFile;
	recordFile.open(fileName);
	if (recordFile.is_open())
	{
		std::cout << "opening " << fileName << " file..." << endl;
		std::string record;
		vector<double> framePoint;
		framePoint.reserve(Dim);
		vector<int> solution;
		while (std::getline(recordFile, record))
		{
			if (frame % timeRange == 0)
			{
				std::istringstream currentRecord(record);
				for (auto i = 0; i < Dim; ++i) {
					double tempNumber = 0;
					currentRecord >> tempNumber;
					framePoint.push_back(tempNumber);
				}
				if (!include_Ball) // remove ball coordinates
				{
					std::rotate(framePoint.begin(), framePoint.begin() + 2, framePoint.end());
					framePoint.pop_back(); framePoint.pop_back();
				}
				// -----------------------
				auto distanceMatrix = getDistanceMatrix(centroids, framePoint);
				auto cost = HungAlg.Solve(distanceMatrix, solution);
				auto refVector = refine_vector(framePoint, solution);

				result_vector.push_back(refVector);
				framePoint.clear();
				solution.clear();
			}
			frame++;
		}
	}
	std::cout << fileName << " extraction complete." << endl;
	recordFile.close();
	return result_vector;
}

static vector<double> getCentroids(std::string filename, int frameRate, int Dim)
{
	auto dataSet = initDataVector(frameRate, filename, Dim, false);
	vector<double> result(dataSet[0].size());

	for (int i = 0; i < dataSet.size(); ++i)
	{
		for (int j = 0; j < result.size(); ++j)
		{
			result[j] += dataSet[i][j];
		}
	}
	for (int i = 0; i < result.size(); ++i)
	{
		result[i] = result[i] / dataSet.size();
	}
	return result;
}
static vector<double> getCentroids(vector<vector<double>>& ds, int frameRate, int Dim)
{

	vector<double> result(ds[0].size());
	for (int i = 0; i < ds.size(); ++i)
	{
		for (int j = 0; j < result.size(); ++j)
		{
			result[j] += ds[i][j];
		}
	}
	for (int i = 0; i < result.size(); ++i)
	{
		result[i] = result[i] / ds.size();
	}
	/*for (int i = 0; i < result.size(); ++i)
	{
		for (int j = 0; j < ds.size(); ++j)
		{
			result[i] += ds[j][i];
		}
	}
	for (int i = 0; i < result.size(); ++i)
	{
		result[i] = result[i] / ds.size();
	}*/
	return result;
}

static void drawVector(vector<double> v, std::string title)
{
	if (title != " ")
	{
		cout << title << endl;
	}
	for (int i = 0; i < v.size(); ++i)
	{
		cout << v[i] << " ";
	}
	cout << endl;
}
static void drawVector(vector<int> v, std::string title)
{
	if (title != " ")
	{
		std::cout << title << endl;
	}
	for (auto i = 0; i < v.size(); ++i)
	{
		std::cout << v[i] << " ";
	}
	std::cout << endl;
}

class FTagRegionInfo
{
public:
	FTagRegionInfo()
		: Tag()
	{}
	FTagRegionInfo(kdtree::vector_t n)
		: Tag(n)
	{}

	kdtree::vector_t Tag;
};
typedef kdtree::KDTree<FTagRegionInfo> KDTree;
// Magic numbers
const double maxR = 30000;
const double minR = -30000;

static std::shared_ptr<KDTree> GetKDtree_with_HA(vector<double> centroids, int timeRange = 75, int Dim = 46, std::string file_name = "", bool include_Ball = true)
{
	std::shared_ptr<KDTree> kd_tree = std::make_shared<KDTree>(Dim - 2);
	ifstream dataSetFile;
	dataSetFile.open(file_name);
	int frame = 0;
	int stepcounter = 0;
	vector<int> assignment;
	HungarianAlgorithm HungAlg;
	// auto centroids = getCentroids(file_name, timeRange, Dim);
	//int howManyFramesinTree = 0;
	// =================== Open data set and build tree
	if (dataSetFile.is_open())
	{
		std::cout << file_name << " data set file opened" << endl;
		kdtree::vector_t hPoint;
		kdtree::vector_t framePoint;
		hPoint.reserve(Dim);
		framePoint.reserve(Dim + 1);
		// Scan data set file and put into kd tree
		// 105000 = 70min
		std::string record;
		while (std::getline(dataSetFile, record))
		{
			frame++;
			if (frame % timeRange == 0)
			{
				std::istringstream currentRecord(record);
				// framePoint.push_back(frame); // TODO: Use it to store frame number of kd-tree node
				for (auto i = 0; i < Dim; ++i) {
					double tempNumber = 0;
					currentRecord >> tempNumber;
					hPoint.push_back(tempNumber);
					framePoint.push_back(tempNumber);
				}
				if (!include_Ball)
				{
					std::rotate(hPoint.begin(), hPoint.begin() + 2, hPoint.end());
					hPoint.pop_back(); hPoint.pop_back();
					// temp
					std::rotate(framePoint.begin(), framePoint.begin() + 2, framePoint.end());
					framePoint.pop_back(); framePoint.pop_back();
					// ----------------					
				}
				try {
					auto costMatrix = getDistanceMatrix(centroids, hPoint);
					auto cost = HungAlg.Solve(costMatrix, assignment);
					auto refined = refine_vector(hPoint, assignment);
					// kd_tree.insert(refined, FTagRegionInfo(framePoint));
					kd_tree->insert(refined, FTagRegionInfo(refined));
					stepcounter++;
					// recordsCollection.push_back(hPoint);
					// howManyFramesinTree++;					
				}
				catch (kdtree::KeyDuplicateException* e) {} // KeyDuplicateException skip adding to tree
				hPoint.clear();
				framePoint.clear();
				assignment.clear();
			}
		}
	}
	dataSetFile.close();
	return kd_tree;
}
static std::shared_ptr<KDTree> GetKDtree(int timeRange = 75, std::string file_name = "", bool include_Ball = true)
{
	// TODO: Refactor Dim calculation depend on include_Ball flag
	int Dim = 46;
	/*if(include_Ball)
	{
		Dim = 46;
	} else
	{
		Dim = 44;
	}*/
	std::shared_ptr<KDTree> kd_tree = std::make_shared<KDTree>(include_Ball ? 46 : 44); // TODO: fix that quick solution
	ifstream dataSetFile;
	dataSetFile.open(file_name);
	int frame = 0;
	int stepcounter = 0;
	// =================== Open data set and build tree
	if (dataSetFile.is_open())
	{
		std::cout << file_name << " data set file opened" << endl;
		kdtree::vector_t hPoint;
		// kdtree::vector_t framePoint;
		hPoint.reserve(Dim);
		// framePoint.reserve(Dim + 1);

		std::string record;
		while (std::getline(dataSetFile, record))
		{
			frame++;
			if (frame % timeRange == 0)
			{
				std::istringstream currentRecord(record);
				// framePoint.push_back(frame); // TODO: Use it to fix frame number of kd-tree node
				for (auto i = 0; i < Dim; ++i) {
					/*if (i != 0 && i != 1)
					{*/
					double tempNumber = 0;
					currentRecord >> tempNumber;
					hPoint.push_back(tempNumber);
					//}
					//framePoint.push_back(tempNumber);
				}

				if (!include_Ball)
				{
					std::rotate(hPoint.begin(), hPoint.begin() + 2, hPoint.end());
					hPoint.pop_back(); hPoint.pop_back();
					// ----------------					
				}
				try {
					kd_tree->insert(hPoint, FTagRegionInfo(hPoint));
					stepcounter++;
					// recordsCollection.push_back(hPoint);
					// howManyFramesinTree++;					
				}
				catch (kdtree::KeyDuplicateException* e)
				{
					cout << "We got duplicate" << endl;
				} // KeyDuplicateException skip adding to tree
				hPoint.clear();
				// framePoint.clear();
			}
		}
	}
	dataSetFile.close();
	return kd_tree;
}

static bool simple_correctness_checker(vector<double> result, vector<double> data_set_point, double range)
{
	for (int i = 0; i < result.size(); ++i)
	{
		auto r = result[i];
		auto t = data_set_point[i];
		auto comp = std::fabs(r - t);

		if (comp > range)
		{
			return false;
		}
	}
	return true;
}
bool isInsideCircle(double x1, double y1, double x2, double y2, double radius)
{
	bool result = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2)) <= radius;
	return result;
}

inline std::pair<kdtree::vector_t, kdtree::vector_t> getTestQuery(kdtree::vector_t queryReference, double range = 1)
{
	kdtree::vector_t lowBound;
	kdtree::vector_t upperBound;
	for (auto i = 0; i < queryReference.size(); ++i) // size = 46
	{
		lowBound.push_back(queryReference[i] - range);
		upperBound.push_back(queryReference[i] + range);
	}
	return std::make_pair(lowBound, upperBound);
}
inline std::pair<kdtree::vector_t, kdtree::vector_t> getTestQueryWithGradient(kdtree::vector_t query_reference)
{
	kdtree::vector_t lowBound;
	kdtree::vector_t upperBound;
	double ballX = query_reference[0];
	double ballY = query_reference[1];
	lowBound.push_back(query_reference[0] - 2);
	upperBound.push_back(query_reference[0] + 2);

	lowBound.push_back(query_reference[1] - 2);
	upperBound.push_back(query_reference[1] + 2);
	for (auto i = 2; i < query_reference.size(); i += 2) // size = 46
	{
		double player_x = query_reference[i];
		double player_y = query_reference[i + 1];
		double distance = sqrt(pow(player_x - ballX, 2) + pow(player_y - ballY, 2));
		// double dyn_range = pow(1.15, distance - 1);
		double dyn_range = pow(1.15, distance + 1);
		// std::cout << cubeRoot << " ";
		lowBound.push_back(query_reference[i] - dyn_range);
		upperBound.push_back(query_reference[i] + dyn_range);

		lowBound.push_back(query_reference[i + 1] - dyn_range);
		upperBound.push_back(query_reference[i + 1] + dyn_range);
	}
	// std::cout << " -------------------- " << std::endl;
	return std::make_pair(lowBound, upperBound);
}
inline std::pair<kdtree::vector_t, kdtree::vector_t> getTestQueryAroundBall(kdtree::vector_t queryReference, double radius, double range = 1)
{
	kdtree::vector_t lowBound;
	kdtree::vector_t upperBound;
	double ballX = queryReference[0];
	double ballY = queryReference[1];

	for (auto i = 0; i < queryReference.size(); ++i) // size = 46
	{
		if (i % 2 == 0)
		{
			if (isInsideCircle(ballX, ballY, queryReference[i], queryReference[i + 1], radius))
			{
				// std::cout << "getTestQueryAroundBall i = " << i << std::endl;
				lowBound.push_back(queryReference[i] - range);
				upperBound.push_back(queryReference[i] + range);

				lowBound.push_back(queryReference[i + 1] - range);
				upperBound.push_back(queryReference[i + 1] + range);
			}
			else
			{
				lowBound.push_back(minR);
				upperBound.push_back(maxR);

				lowBound.push_back(minR);
				upperBound.push_back(maxR);
			}
		}
	}
	return std::make_pair(lowBound, upperBound);
}

static int getBestResult(vector<double> origin, vector<vector<double>> resultSet)
{
	double sum = 0;
	vector<double> sums(resultSet.size());
	for (int j = 0; j < resultSet.size(); ++j)
	{
		for (int i = 2; i < origin.size(); i += 2)
		{
			auto x1 = origin[i];
			auto y1 = origin[i + 1];

			auto x2 = resultSet[j][i];
			auto y2 = resultSet[j][i + 1];

			sum += sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
		}
		sums[j] = sum;
		sum = 0;
	}
	std::vector<double>::iterator result = std::min_element(sums.begin(), sums.end());
	return std::distance(sums.begin(), result);
	/*cout << "min element at: " << std::distance(sums.begin(), result) << endl;
	cout << "min is " << sums[std::distance(sums.begin(), result)] << endl;*/
}
static int getBestResult(vector<double> origin, vector<FTagRegionInfo> resultSet)
{
	double sum = 0;
	vector<double> sums(resultSet.size());
	for (int j = 0; j < resultSet.size(); ++j)
	{
		for (int i = 2; i < origin.size(); i += 2)
		{
			auto x1 = origin[i];
			auto y1 = origin[i + 1];
						
			auto x2 = resultSet[j].Tag[i];
			auto y2 = resultSet[j].Tag[i + 1];

			sum += sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
		}
		sums[j] = sum;
		sum = 0;
	}
	std::vector<double>::iterator result = std::min_element(sums.begin(), sums.end());
	return std::distance(sums.begin(), result);
}