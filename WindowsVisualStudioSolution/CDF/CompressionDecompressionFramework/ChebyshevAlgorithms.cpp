#include "ChebyshevAlgorithms.h"

// Type alias for cache key, where first is the degree n and second is the value x
using CacheKey = std::pair<int, double>;

// Custom hash function for CacheKey
struct CacheKeyHash {
	std::size_t operator()(const CacheKey& key) const {
		return std::hash<int>()(key.first) ^ std::hash<double>()(key.second);
	}
};

std::unordered_map<CacheKey, double, CacheKeyHash> cache;

double ChebyshevAlgorithms::chebyshevT(int n, double x) {
	CacheKey key = {n, x};
	
	// Check if value is in cache
	auto iter = cache.find(key);
	if (iter != cache.end()) {
		return iter->second;
	}

	// Calculate Chebyshev polynomial
	double result;
	if (n == 0) {
		result = 1.0;
	} else if (n == 1) {
		result = x;
	} else {
		double T_prev = 1.0;
		double T_curr = x;
		for (int i = 2; i <= n; ++i) {
			double T_next = 2 * x * T_curr - T_prev;
			T_prev = T_curr;
			T_curr = T_next;
		}
		result = T_curr;
	}

	// Store result in cache and return
	cache[key] = result;
	return result;
}

