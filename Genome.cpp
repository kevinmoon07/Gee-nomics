#include "provided.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

class GenomeImpl
{
public:
	GenomeImpl(const string& nm, const string& sequence);
	static bool load(istream& genomeSource, vector<Genome>& genomes);
	int length() const;
	string name() const;
	bool extract(int position, int length, string& fragment) const;
private:
	vector<char> m_name;
	vector<char> m_sequence;
};

GenomeImpl::GenomeImpl(const string& nm, const string& sequence)
{
	for (int x = 0; x < nm.length(); x++)
	{
		m_name.push_back(nm[x]);
	}

	for (int x = 0; x < sequence.length(); x++)
	{
		m_sequence.push_back(sequence[x]);
	}
}

bool GenomeImpl::load(istream& genomeSource, vector<Genome>& genomes)
{
	genomes.clear();

	string sequence = "";
	string name = "";
	char ch;
	genomeSource.get(ch);
	bool endNewLine = true;

	//Return false if the file is not properly formatted
	if (!genomeSource)
	{
		return false;
	}

	bool charAfterName = false;

	//Return false if the first character is not the beginning of a name
	//The file must start with a name line
	if (ch == '>')
	{

		//Store the first genome's name
		while (genomeSource.get(ch) && ch != '\n')
		{
			if (ch != ' ')
			{
				charAfterName = true;
			}

			name += ch;
		}
	}
	else
	{
		return false;
	}

	//If the name does not have any characters in it, return false
	if (!charAfterName)
	{
		return false;
	}

	bool gotBaseAfterName = false;
	genomeSource.get(ch);
	//Loop through each character in the file
	do
	{
		//If the character is the start of a new line
		if (ch == '\n')
		{
			//If it's an empty line, return false
			if (endNewLine)
			{
				return false;
			}

			endNewLine = true;
		}
		else if (ch == ' ')
		{
			endNewLine = false;
		}
		//If the character is the start of a name
		else if (ch == '>' && endNewLine)
		{
			//Return false if there were no base lines after a name line
			if (!gotBaseAfterName)
			{
				return false;
			}

			gotBaseAfterName = false;
			charAfterName = false;

			//Create and store the previous genome in the vector
			Genome x(name, sequence);
			genomes.push_back(x);
			name = "";
			sequence = "";

			//Store the new name
			while (genomeSource.get(ch) && ch != '\n')
			{
				if (ch != ' ')
				{
					charAfterName = true;
				}

				name += ch;
			}

			if (!charAfterName)
			{
				return false;
			}
			endNewLine = true;
		}

		//If the character is a valid character
		else if (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'N'
			|| ch == 'a' || ch == 'c' || ch == 'g' || ch == 't' || ch == 'n')
		{
			//Change the character to uppercase if needed then concatonate it to the sequence
			if (ch == 'a' || ch == 'c' || ch == 'g' || ch == 't' || ch == 'n')
			{
				ch -= 32;
			}
			sequence += ch;
			endNewLine = false;
			gotBaseAfterName = true;
		}

		//If the character is not a valid character, return false
		else
		{
			return false;
		}

	} while (genomeSource.get(ch));


	//Add in the last genome if it had a sequence
	if (gotBaseAfterName)
	{
		Genome x(name, sequence);
		genomes.push_back(x);
		return true;
	}
	else
	{
		return false;
	}
}

int GenomeImpl::length() const
{
	return m_sequence.size();
}

string GenomeImpl::name() const
{
	string name = "";
	for (int x = 0; x < m_name.size(); x++)
	{
		name += m_name[x];
	}

	return name;
}

bool GenomeImpl::extract(int position, int length, string& fragment) const
{
	//Return false if user is trying to return a subset that goes past the length of the genome
	if (position + length > m_sequence.size())
	{
		return false;
	}

	//Return false if the position is past the fragment length
	if (position >= m_sequence.size())
	{
		return false;
	}

	string fragmentTemp = "";

	//Iterate through the genome at the specified starting location for the specified length
	for (int posTemp = position; posTemp < position + length; posTemp++)
	{
		fragmentTemp += m_sequence[posTemp];
	}

	//Change the value of fragment
	fragment = fragmentTemp;

	return true;
}

//******************** Genome functions ************************************

// These functions simply delegate to GenomeImpl's functions.
// You probably don't want to change any of this code.

Genome::Genome(const string& nm, const string& sequence)
{
	m_impl = new GenomeImpl(nm, sequence);
}

Genome::~Genome()
{
	delete m_impl;
}

Genome::Genome(const Genome& other)
{
	m_impl = new GenomeImpl(*other.m_impl);
}

Genome& Genome::operator=(const Genome& rhs)
{
	GenomeImpl* newImpl = new GenomeImpl(*rhs.m_impl);
	delete m_impl;
	m_impl = newImpl;
	return *this;
}

bool Genome::load(istream& genomeSource, vector<Genome>& genomes)
{
	return GenomeImpl::load(genomeSource, genomes);
}

int Genome::length() const
{
	return m_impl->length();
}

string Genome::name() const
{
	return m_impl->name();
}

bool Genome::extract(int position, int length, string& fragment) const
{
	return m_impl->extract(position, length, fragment);
}
