#include "provided.h"
#include "Trie.h"
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>

using namespace std;

class GenomeMatcherImpl
{
public:
	GenomeMatcherImpl(int minSearchLength);
	void addGenome(const Genome& genome);
	int minimumSearchLength() const;
	bool findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
	bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const;
private:
	int m_minSearchLength;
	vector<Genome> m_genomes;
//First int is the position of the genome in question in  the vector starting at 0
	//Second int is the position of the start of the sequence that ends the Trie
	Trie<pair<int, int>> m_trie; 
};

GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength)
{
	m_minSearchLength = minSearchLength;
}

void GenomeMatcherImpl::addGenome(const Genome& genome)
{
	m_genomes.push_back(genome);
	string shortSeq = ""; //Create a blank string to store the subsequence
	pair<int, int> specimen;
	specimen.first = m_genomes.size() - 1;

	for (int pos = 0; pos <= genome.length() - m_minSearchLength; pos++)
	{
		genome.extract(pos, m_minSearchLength, shortSeq); //Pass in the subsequence
		specimen.second = pos;

		m_trie.insert(shortSeq, specimen); //Insert the new subsequence and Value into the Trie
		shortSeq = ""; //Reset the subsequence
	}
}

int GenomeMatcherImpl::minimumSearchLength() const
{
	return m_minSearchLength;
}

bool GenomeMatcherImpl::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
	if (fragment.length() < minimumLength || minimumLength < m_minSearchLength)
	{
		return false;
	}

	//Make sure key is valid
	for (int checkCase = 0; checkCase < fragment.length(); checkCase++)
	{
		if (!isalpha(fragment[checkCase]))
		{
			return false;
		}
		else if (!isupper(fragment[checkCase]))
		{
			return false;
		}
	}

	vector<pair<int, int>> trieFind;

	int extractionLength = fragment.length() - m_minSearchLength;
	//If extraction Length is 0, that means that fragment length and minSearchLength are equal, so no bases need to
	//be extracted
	if (extractionLength == 0)
	{
		trieFind = m_trie.find(fragment, exactMatchOnly);
		if (trieFind.empty())
		{
			return false;
		}

		//Sort by index first in increasing order, then by pos in increasing order
		sort(trieFind.begin(), trieFind.end());

		vector<pair<int, int>>::iterator it;
		it = trieFind.begin();
		vector<pair<int, int>>::iterator it2;

		while (it2 != trieFind.end())
		{
			it2 = it + 1;

			//Remove the duplicates and higher position entries
			if (it->first == it2->first)
			{
				if (it->second <= it2->second)
				{
					trieFind.erase(it2);
				}
			}

			//Update the first iterator a position after all the duplicates of a particular genome have been deleted
			else
			{
				it++;
			}

			//Set it2 to what it will be evaluating in the next iteration to see if loop should end
			it2 = it + 1;
		}

		//Update the matches
		for (int x = 0; x < trieFind.size(); x++)
		{
			DNAMatch g;
			g.genomeName = m_genomes[trieFind[x].first].name();
			g.position = trieFind[x].second;
			g.length = fragment.length();

			matches.push_back(g);
		}
		return true;
	}

	//Extraction length is greater than 0
	else
	{
		string fragmentTemp = fragment.substr(0, m_minSearchLength);

		//Return potential locations
		trieFind = m_trie.find(fragmentTemp, exactMatchOnly);

		//Return false if no potential locations
		if (trieFind.empty())
		{
			return false;
		}

		sort(trieFind.begin(), trieFind.end());

		vector<pair<int, int>>::iterator it;
		it = trieFind.begin();
		vector<pair<int, int>>::iterator it2;
		it2 = it + 1;

		//Delete all duplicating genomes
		while (it2 != trieFind.end())
		{

			it2 = it + 1;

			//Remove the duplicates
			if (it->first == it2->first)
			{
				if (it->second == it2->second)
				{
					trieFind.erase(it2);
				}
				else
				{
					it++;
				}
			}

			//Update the first iterator a position after all the duplicates of a particular genome have been deleted
			else
			{
				it++;
			}

			//Set it2 to what it will be evaluating in the next iteration to see if loop should end
			it2 = it + 1;
		}

		string fragmentExtracted = "";

		vector<Genome>::const_iterator itG;

		//Go through each trieFind element
		for (int x = 0; x < trieFind.size(); x++)
		{

			////Find the genome in the library that has the same name
			//itG = m_genomes.begin();

			//while (trieFind[x].first != itG->name())
			//{
			//	itG++;
			//	if (itG == m_genomes.end())
			//	{
			//		return false;
			//	}
			//}

			//Extract the entire sequence
			m_genomes[trieFind[x].first].extract(trieFind[x].second, fragment.length(), fragmentExtracted);

			//Compare the fragment and extracted fragment
			bool matchMaker = exactMatchOnly;
			int y = 0;
			for (y = 0; y < fragment.length(); y++)
			{
				//If the two fragments are not exact matches, break out of the loop
				if (matchMaker)
				{
					if (fragment[y] != fragmentExtracted[y])
					{
						break;
					}
				}

				//If the two characters are equal and matchMaker is false, keep it as such. Otherwise, switch
				else
				{
					if (fragment[y] != fragmentExtracted[y])
					{
						matchMaker = true;
					}
				}
			}

			int length = y;

			//If valid, add the DNA Match to matches
			if (length >= minimumLength)
			{
				DNAMatch g;
				g.genomeName = m_genomes[trieFind[x].first].name();
				g.position = trieFind[x].second;
				g.length = length;

				matches.push_back(g);
			}

			fragmentExtracted = "";
		}
		return true;
	}
}

bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
	if (fragmentMatchLength < m_minSearchLength)
	{
		return false;
	}

	vector<int> sequencePos;
	string fragment = "";
	vector<DNAMatch> matches;
	int S = query.length() / fragmentMatchLength;

	map<string, int> foundMatchCount;
	map<string, int>::iterator it;

	for (int x = 0; x < S; x++)
	{
		sequencePos.push_back(x*fragmentMatchLength);
		query.extract(sequencePos[x], fragmentMatchLength, fragment);
		findGenomesWithThisDNA(fragment, fragmentMatchLength, exactMatchOnly, matches);

		for (int y = 0; y < matches.size(); y++)
		{
			//If the genome doesn't exist in the map, create it
			if (foundMatchCount.find(matches[y].genomeName) == foundMatchCount.end())
			{
				foundMatchCount[matches[y].genomeName] = 1;
			}

			//If it does exist, increment the match count
			else
			{
				it = foundMatchCount.find(matches[y].genomeName);
				it->second += 1;
			}
		}
	}

	double percentage;
	results.clear();

	for (it = foundMatchCount.begin(); it != foundMatchCount.end(); it++)
	{
		double dec = static_cast<double>(it->second) / S;
		percentage = dec * 100;
		if (percentage >= matchPercentThreshold)
		{
			GenomeMatch g;
			g.genomeName = it->first;
			g.percentMatch = percentage;
			results.push_back(g);
		}
	}

	if (results.empty())
	{
		return false;
	}
	else
	{
		//Sort results Selection Sort
		for (int x = 0; x < results.size(); x++)
		{
			int minIndex = x;
			for (int y = x + 1; y < results.size(); y++)
			{
				if (results[y].percentMatch > results[minIndex].percentMatch)
				{
					minIndex = y;
				}
				else if (results[y].percentMatch == results[minIndex].percentMatch)
				{
					if (results[y].genomeName < results[minIndex].genomeName)
					{
						minIndex = y;
					}
				}
			}
			swap(results[x], results[minIndex]);
		}
		return true;
	}
}

//******************** GenomeMatcher functions ********************************

// These functions simply delegate to GenomeMatcherImpl's functions.
// You probably don't want to change any of this code.

GenomeMatcher::GenomeMatcher(int minSearchLength)
{
	m_impl = new GenomeMatcherImpl(minSearchLength);
}

GenomeMatcher::~GenomeMatcher()
{
	delete m_impl;
}

void GenomeMatcher::addGenome(const Genome& genome)
{
	m_impl->addGenome(genome);
}

int GenomeMatcher::minimumSearchLength() const
{
	return m_impl->minimumSearchLength();
}

bool GenomeMatcher::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
	return m_impl->findGenomesWithThisDNA(fragment, minimumLength, exactMatchOnly, matches);
}

bool GenomeMatcher::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
	return m_impl->findRelatedGenomes(query, fragmentMatchLength, exactMatchOnly, matchPercentThreshold, results);
}
