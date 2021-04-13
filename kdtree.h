#pragma once

#include <vector>
#include <limits.h>
#include <string>
#include <cmath>
namespace kdtree
{
	//typedef float value_type;
	typedef double value_type;
	typedef std::vector<value_type> vector_t;

#define NegativeInfinity FLT_MIN
#define PositiveInfinity FLT_MAX

#define REMOVE_HIGHEST  1
#define REMOVE_LOWEST  2

	inline value_type sqrt(value_type v)
	{
		
		return std::sqrt(v);
	}
	inline value_type Max(value_type x, value_type y)
	{
		return x > y ? x : y;
	}
	inline value_type Min(value_type x, value_type y)
	{
		return x < y ? x : y;
	}


	class HPoint
	{
	public:
		vector_t coord;
		HPoint(int n) : coord(n)
		{}

		HPoint(const vector_t& x)
			: coord(x)
		{}

		HPoint clone()
		{
			return *this;
		}

		bool equals(const HPoint& p) const
		{
			if (p.coord.size() == coord.size())
			{
				for (int i = 0; i < p.coord.size(); i++)
				{
					if (p.coord[i] != coord[i])
					{
						return false;
					}
				}
				return true;
			}
			return false;
		}
		static value_type sqrdist(const HPoint& x, const HPoint& y)
		{
			value_type dist = 0;
			for (int i = 0; i < x.coord.size(); i++)
			{
				value_type diff = (x.coord[i] - y.coord[i]);
				dist += diff * diff;
			}
			return dist;
		}
		static value_type eucdist(const HPoint& x, const HPoint& y)
		{
			return sqrt(sqrdist(x, y));
		}
	};
	class HRect
	{
	public:
		HPoint min;
		HPoint max;

		HRect() : min(0), max(0)
		{}

		HRect clone() const
		{
			return *this;
		}
		HPoint closest(const HPoint& t)
		{
			HPoint p(t.coord.size());
			for (int i = 0; i < t.coord.size(); ++i)
			{
				if (t.coord[i] <= min.coord[i])
				{
					p.coord[i] = min.coord[i];
				}
				else if (t.coord[i] >= max.coord[i])
				{
					p.coord[i] = max.coord[i];
				}
				else
				{
					p.coord[i] = t.coord[i];
				}
			}
			return p;
		}
		static HRect infiniteHRect(int d)
		{

			HPoint vmin(d);
			HPoint vmax(d);

			for (int i = 0; i < d; ++i)
			{
				vmin.coord[i] = NegativeInfinity;
				vmax.coord[i] = PositiveInfinity;
			}

			return HRect(vmin, vmax);
		}
		HRect intersection(const HRect& r)
		{
			HPoint newmin(min.coord.size());
			HPoint newmax(min.coord.size());

			for (int i = 0; i < min.coord.size(); ++i)
			{
				newmin.coord[i] = Max(min.coord[i], r.min.coord[i]);
				newmax.coord[i] = Min(max.coord[i], r.max.coord[i]);
				if (newmin.coord[i] >= newmax.coord[i]) return HRect();
			}

			return HRect(newmin, newmax);
		}
		value_type area()
		{

			value_type a = 1;

			for (int i = 0; i < min.coord.size(); ++i)
			{
				a *= (max.coord[i] - min.coord[i]);
			}

			return a;
		}
	protected:
		HRect(int ndims) : min(ndims), max(ndims)
		{}
		HRect(const HPoint& vmin, const HPoint& vmax)
			: min(vmin), max(vmax)
		{}
	};

	class KeyDuplicateException
	{};
	class KeySizeException
	{};
	class KeyMissingException
	{};
	class ArgumentException
	{
	public:
		ArgumentException(const char* msg) {}
	};


	template<class T>
	class PriorityQueue
	{
	private:
		std::vector<T*> data;
		std::vector<value_type> value;
		int count;
		int capacity;
	public:
		PriorityQueue()
		{
			init(20, PositiveInfinity);
		}
		PriorityQueue(int capacity)
		{
			init(capacity, PositiveInfinity);
		}
		PriorityQueue(int capacity, value_type maxPriority)
		{
			init(capacity, maxPriority);
		}
		void init(int size, value_type maxPriority)
		{
			capacity = size;
			for (int i = 0; i < capacity + 1; i++)
			{
				data.push_back(NULL);
				value.push_back(0);
			}
			value[0] = maxPriority;
			data[0] = NULL;
			count = 0;
		}
		void add(T* element, value_type priority)
		{
			if (count++ >= capacity)
			{
				expandCapacity();
			}
			/* put this as the last element */
			value[count] = priority;
			data[count] = element;
			bubbleUp(count);
		}
		T* remove()
		{
			if (count == 0)
				return NULL;
			T* element = data[1];
			/* swap the last element into the first */
			data[1] = data[count];
			value[1] = value[count];
			/* let the GC clean up */
			data[count] = NULL;
			value[count] = 0;
			count--;
			bubbleDown(1);
			return element;
		}
		T* front()
		{
			return data[1];
		}
		value_type getMaxPriority()
		{
			return value[1];
		}

		void clear()
		{
			for (int i = 1; i < count; i++)
			{
				data[i] = NULL; /* help gc */
			}
			count = 0;
		}
		int length()
		{
			return count;
		}
	private:
		void bubbleDown(int pos)
		{
			T* element = data[pos];
			value_type priority = value[pos];
			int child;
			/* hole is position '1' */
			for (; pos * 2 <= count; pos = child)
			{
				child = pos * 2;
				/* if 'child' equals 'count' then there
				   is only one leaf for this parent */
				if (child != count)

					/* left_child > right_child */
					if (value[child] < value[child + 1])
						child++; /* choose the biggest child */
				/* percolate down the data at 'pos', one level
				   i.e biggest child becomes the parent */
				if (priority < value[child])
				{
					value[pos] = value[child];
					data[pos] = data[child];
				}
				else
				{
					break;
				}
			}
			value[pos] = priority;
			data[pos] = element;
		}
		void bubbleUp(int pos)
		{
			T* element = data[pos];
			value_type priority = value[pos];
			/* when the parent is not less than the child, end*/
			while (value[pos / 2] < priority)
			{
				/* overwrite the child with the parent */
				value[pos] = value[pos / 2];
				data[pos] = data[pos / 2];
				pos /= 2;
			}
			value[pos] = priority;
			data[pos] = element;
		}
		void expandCapacity()
		{
			capacity = count * 2;
			for (int i = 0; i < capacity + 1; i++)
			{
				data.push_back(NULL);
				value.push_back(0);
			}
		}

	};

	template<class T>
	class NearestNeighborList
	{
	public:
		PriorityQueue<T> m_Queue;
		int m_Capacity;
		NearestNeighborList(int capacity)
			: m_Capacity(capacity), m_Queue(capacity, PositiveInfinity)
		{
		}
		value_type getMaxPriority()
		{
			if (m_Queue.length() == 0)
			{
				return PositiveInfinity;
			}
			return m_Queue.getMaxPriority();
		}
		bool insert(T* _object, value_type priority)
		{
			if (m_Queue.length() < m_Capacity)
			{
				// capacity not reached
				m_Queue.add(_object, priority);
				return true;
			}
			if (priority > m_Queue.getMaxPriority())
			{
				// do not insert - all elements in queue have lower priority
				return false;
			}
			// remove object with highest priority
			m_Queue.remove();
			// add new object
			m_Queue.add(_object, priority);
			return true;
		}
		bool isCapacityReached()
		{
			return m_Queue.length() >= m_Capacity;
		}
		T* getHighest()
		{
			return m_Queue.front();
		}
		bool isEmpty()
		{
			return m_Queue.length() == 0;
		}
		int getSize()
		{
			return m_Queue.length();
		}
		T* removeHighest()
		{
			// remove object with highest priority
			return m_Queue.remove();
		}
	};


	// 附带信息
	class DefaultKDTag
	{
	public:
		DefaultKDTag(std::string val)
		{
			Val = val;
		}
		DefaultKDTag()
		{
			
		}
		std::string Val = "this is 00";
	};
	template<class KDTag>
	class TkdNode
	{
	protected:
		HPoint k;
		TkdNode* left, *right;

		TkdNode(const HPoint& key, const KDTag& val)
			: k(key), left(NULL), right(NULL), v(val), deleted(false)
		{}

	public:
		KDTag v;
		bool deleted;

		static void walk(TkdNode* root, std::vector<TkdNode*>& walker)
		{
			if (root == NULL)
				return;

			walker.push_back(root);
			walk(root->left, walker);
			walk(root->right, walker);
		}
		static void recursive_delete(TkdNode* root)
		{
			std::vector<TkdNode*> walker;
			TkdNode::walk(root, walker);
			for (int i = 0; i < walker.size(); i++)
			{
				delete walker[i];
			}
			walker.clear();
		}

		static TkdNode* ins(const HPoint& key, const KDTag& val, TkdNode* t, int lev, int K)
		{
			if (t == NULL)
			{
				t = new TkdNode(key, val);
			}

			else if (key.equals(t->k))
			{

				// "re-insert"
				if (t->deleted)
				{
					t->deleted = false;
					t->v = val;
				}

				else
				{
					throw (new KeyDuplicateException());
				}
			}

			else if (key.coord[lev] > t->k.coord[lev])
			{
				t->right = ins(key, val, t->right, (lev + 1) % K, K);
			}
			else
			{
				t->left = ins(key, val, t->left, (lev + 1) % K, K);
			}

			return t;
		}
		static TkdNode* srch(const HPoint& key, TkdNode* t, int K)
		{
			for (int lev = 0; t != NULL; lev = (lev + 1) % K)
			{

				if (!t->deleted && key.equals(t->k))
				{
					return t;
				}
				else if (key.coord[lev] > t->k.coord[lev])
				{
					t = t->right;
				}
				else
				{
					t = t->left;
				}
			}

			return NULL;
		}
		static TkdNode* Delete(const HPoint& key, TkdNode* t, int lev, int K, bool& deleted) {
			if (t == NULL) return NULL;

			if (!t->deleted && key.equals(t->k))
			{
				t->deleted = true;
				deleted = true;
			}
			else if (key.coord[lev] > t->k.coord[lev])
			{
				t->right = Delete(key, t->right, (lev + 1) % K, K, deleted);
			}
			else
			{
				t->left = Delete(key, t->left, (lev + 1) % K, K, deleted);
			}

			if (!t->deleted || t->left != NULL || t->right != NULL)
			{
				return t;
			}
			else
			{
				return NULL;
			}
		}
		// range search itself
		static void rsearch(const HPoint& lowk, const HPoint& uppk, TkdNode* t, int lev,
			int K, std::vector<TkdNode*>& v)
		{

			if (t == NULL) return;
			if (lowk.coord[lev] <= t->k.coord[lev])
			{
				rsearch(lowk, uppk, t->left, (lev + 1) % K, K, v);
			}
			int j;
			for (j = 0; j < K && lowk.coord[j] <= t->k.coord[j] &&
				uppk.coord[j] >= t->k.coord[j]; j++)
				;
			if (j == K && !t->deleted) v.push_back(t);
			if (uppk.coord[lev] > t->k.coord[lev])
			{
				rsearch(lowk, uppk, t->right, (lev + 1) % K, K, v);
			}
		}
		static void nnbr(TkdNode* kd, const HPoint& target, HRect& hr,
			value_type max_dist_sqd, int lev, int K,
			NearestNeighborList<TkdNode>& nnl)
		{

			// 1. if kd is empty then set dist-sqd to infinity and exit.
			if (kd == NULL)
			{
				return;
			}

			// 2. s := split field of kd
			int s = lev % K;

			// 3. pivot := dom-elt field of kd
			const HPoint& pivot = kd->k;
			value_type pivot_to_target = HPoint::sqrdist(pivot, target);

			// 4. Cut hr into to sub-hyperrectangles left-hr and right-hr.
			//    The cut plane is through pivot and perpendicular to the s
			//    dimension.
			HRect& left_hr = hr; // optimize by not cloning
			HRect right_hr = hr.clone();
			left_hr.max.coord[s] = pivot.coord[s];
			right_hr.min.coord[s] = pivot.coord[s];

			// 5. target-in-left := target_s <= pivot_s
			bool target_in_left = target.coord[s] < pivot.coord[s];

			TkdNode* nearer_kd;
			HRect& nearer_hr = right_hr;
			TkdNode* further_kd;
			HRect& further_hr = left_hr;

			// 6. if target-in-left then
			//    6.1. nearer-kd := left field of kd and nearer-hr := left-hr
			//    6.2. further-kd := right field of kd and further-hr := right-hr
			if (target_in_left)
			{
				nearer_kd = kd->left;
				nearer_hr = left_hr;
				further_kd = kd->right;
				further_hr = right_hr;
			}
			//
			// 7. if not target-in-left then
			//    7.1. nearer-kd := right field of kd and nearer-hr := right-hr
			//    7.2. further-kd := left field of kd and further-hr := left-hr
			else
			{
				nearer_kd = kd->right;
				nearer_hr = right_hr;
				further_kd = kd->left;
				further_hr = left_hr;
			}

			// 8. Recursively call Nearest Neighbor with paramters
			//    (nearer-kd, target, nearer-hr, max-dist-sqd), storing the
			//    results in nearest and dist-sqd
			nnbr(nearer_kd, target, nearer_hr, max_dist_sqd, lev + 1, K, nnl);

			TkdNode* nearest = nnl.getHighest();
			value_type dist_sqd;

			if (!nnl.isCapacityReached())
			{
				dist_sqd = PositiveInfinity;
			}
			else
			{
				dist_sqd = nnl.getMaxPriority();
			}

			// 9. max-dist-sqd := minimum of max-dist-sqd and dist-sqd
			max_dist_sqd = Min(max_dist_sqd, dist_sqd);

			// 10. A nearer point could only lie in further-kd if there were some
			//     part of further-hr within distance sqrt(max-dist-sqd) of
			//     target.  If this is the case then
			HPoint closest = further_hr.closest(target);
			if (HPoint::eucdist(closest, target) < sqrt(max_dist_sqd))
			{

				// 10.1 if (pivot-target)^2 < dist-sqd then
				if (pivot_to_target < dist_sqd)
				{

					// 10.1.1 nearest := (pivot, range-elt field of kd)
					nearest = kd;

					// 10.1.2 dist-sqd = (pivot-target)^2
					dist_sqd = pivot_to_target;

					// add to nnl
					if (!kd->deleted)
					{
						nnl.insert(kd, dist_sqd);
					}

					// 10.1.3 max-dist-sqd = dist-sqd
					// max_dist_sqd = dist_sqd;
					if (nnl.isCapacityReached())
					{
						max_dist_sqd = nnl.getMaxPriority();
					}
					else
					{
						max_dist_sqd = PositiveInfinity;
					}
				}

				// 10.2 Recursively call Nearest Neighbor with parameters
				//      (further-kd, target, further-hr, max-dist_sqd),
				//      storing results in temp-nearest and temp-dist-sqd
				nnbr(further_kd, target, further_hr, max_dist_sqd, lev + 1, K, nnl);
				TkdNode* temp_nearest = nnl.getHighest();
				value_type temp_dist_sqd = nnl.getMaxPriority();

				// 10.3 If tmp-dist-sqd < dist-sqd then
				if (temp_dist_sqd < dist_sqd)
				{

					// 10.3.1 nearest := temp_nearest and dist_sqd := temp_dist_sqd
					nearest = temp_nearest;
					dist_sqd = temp_dist_sqd;
				}
			}

			// SDL: otherwise, current point is nearest
			else if (pivot_to_target < max_dist_sqd)
			{
				nearest = kd;
				dist_sqd = pivot_to_target;
			}
		}
	};

	template<class KDTag>
	class KDTree
	{
		typedef TkdNode<KDTag> KDNode;
	protected:
		int m_K;
		KDNode* m_root;
		int m_count;

	public:
		KDTree(int k)
		{
			m_K = k;
			m_root = NULL;
			m_count = 0;
		}
		~KDTree()
		{
			Reset();
		}
		bool IsValid()
		{
			return m_count > 0;
		}
		void Reset()
		{
			// 根据m_root，销毁所有的KDNode
			KDNode::recursive_delete(m_root);
		}
		void insert(vector_t key, KDTag value)
		{

			if (key.size() != m_K)
			{
				throw new KeySizeException();
			}

			try
			{
				m_root = KDNode::ins(HPoint(key), value, m_root, 0, m_K);
			}

			catch (KeyDuplicateException e)
			{
				throw e;
			}

			m_count++;
		}
		KDTag search(vector_t key)
		{

			if (key.size() != m_K)
			{
				throw new KeySizeException();
			}

			KDNode* kd = KDNode::srch(HPoint(key), m_root, m_K);

			return (kd == NULL ? KDTag() : kd->v);
		}
		void Delete(vector_t key)
		{

			if (key.size() != m_K)
			{
				throw new KeySizeException();
			}

			else
			{
				bool deleted = false;
				m_root = KDNode::Delete(HPoint(key), m_root, 0, m_K, deleted);
				if (deleted == false) {
					throw new KeyMissingException();
				}
				m_count--;
			}
		}
		KDTag nearest(const vector_t& key) const
		{
			std::vector<KDTag> nbrs = nearest(key, 1);
			return nbrs[0];
		}
		std::vector<KDTag> nearest(const vector_t& key, int n) const
		{
			if (n < 0 || n > m_count)
			{
				throw new ArgumentException("Number of neighbors cannot be negative or greater than number of nodes");
			}

			if (key.size() != m_K)
			{
				throw new KeySizeException();
			}

			std::vector<KDTag> nbrs(n);
			NearestNeighborList<KDNode> nnl(n);

			// initial call is with infinite hyper-rectangle and max distance
			HRect hr = HRect::infiniteHRect(key.size());
			value_type max_dist_sqd = PositiveInfinity;
			HPoint keyp(key);

			KDNode::nnbr(m_root, keyp, hr, max_dist_sqd, 0, m_K, nnl);

			for (int i = 0; i < n; ++i)
			{
				KDNode* kd = nnl.removeHighest();
				nbrs[n - i - 1] = kd->v;
			}

			return nbrs;
		}
		std::vector<KDTag> range(const vector_t& lowk, const vector_t& uppk)
		{

			if (lowk.size() != uppk.size())
			{
				throw new KeySizeException();
			}

			else if (lowk.size() != m_K)
			{
				throw new KeySizeException();
			}

			else
			{
				std::vector<KDNode*> v;
				KDNode::rsearch(HPoint(lowk), HPoint(uppk),
					m_root, 0, m_K, v);
				std::vector<KDTag> o(v.size());
				for (int i = 0; i < v.size(); ++i)
				{
					KDNode* n = v[i];
					o[i] = n->v;
				}
				return o;
			}
		}
	};

	class TestKDTree
	{
	public:
		void Test()
		{
			KDTree<DefaultKDTag> kdtreeYiZhan(2);
			vector_t key1 = { 1, 1 };
			std::string val1;
			for (int i = 0; i < key1.size(); ++i)
			{
				val1.append(std::to_string(key1[i]));
			}
			kdtreeYiZhan.insert(key1, DefaultKDTag(val1));
			//-------------
			vector_t key2 = { 2, 2 };
			std::string val2;
			for (int i = 0; i < key2.size(); ++i)
			{
				val2.append(std::to_string(key2[i]));
			}
			kdtreeYiZhan.insert(key2, DefaultKDTag(val2));
			DefaultKDTag tag = kdtreeYiZhan.nearest({ 3, 3 });
		}
	};
}


