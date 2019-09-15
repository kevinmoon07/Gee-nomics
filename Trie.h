#ifndef TRIE_INCLUDED
#define TRIE_INCLUDED

#include <string>
#include <vector>


using namespace std;

template<typename ValueType>
class Trie
{
public:
	Trie(); //Constructs a new Node
	~Trie(); //Deletes all existing nodes
	void reset(); //Deletes all existing nodes and creates a new node
	void insert(const std::string& key, const ValueType& value);
	std::vector<ValueType> find(const std::string& key, bool exactMatchOnly) const;

	// C++11 syntax for preventing copying and assignment
	Trie(const Trie&) = delete;
	Trie& operator=(const Trie&) = delete;
private:
	vector<ValueType> findHelper(const std::string& key, bool exactMatchOnly, const Trie* pos) const;
	vector<Trie*> children;
	vector<ValueType> values;
	Trie* m_root;

};

template<typename ValueType>
Trie<ValueType>::Trie()
{
	children = vector<Trie*>(26); //Allocate a vector of Tries
	m_root = this;
}

template<typename ValueType>
Trie<ValueType>::~Trie()
{
	for (int count = 0; count < 26; count++) //Loop through each child vector and pass through the subsequent node
	{
		if (this->children[count] != nullptr) //If the vector is not null, pass through the subsequent node
		{
			delete this->children[count];
		}
	}
}

template<typename ValueType>
void Trie<ValueType>::reset()
{
	delete m_root;
	m_root = new Trie;
}

template<typename ValueType>
void Trie<ValueType>::insert(const string& key, const ValueType& value)
{
	Trie* pos = this;
	for (int x = 0; x < key.length(); x++) //Loop through the word and...
	{
		if (pos->children[key[x] - 'A'] == nullptr)
			//If the character has not been initialized, initialize it with a new Trie
		{
			pos->children[key[x] - 'A'] = new Trie;
		}			
		pos = pos->children[key[x] - 'A']; //Update the position pointer to the next Trie
	}

	pos->values.push_back(value); //At the end position, insert value into the values vector
}

template<typename ValueType>
vector<ValueType> Trie<ValueType>::findHelper(const std::string& key, bool exactMatchOnly, const Trie* pos) const 
{
	vector<ValueType> totValue;
	string nextKey = key.substr(1, key.length() - 1);
	if (exactMatchOnly)
	{
		for (int x = 0; x < key.length(); x++) //Loop through each letter in the key and update the position
		{
			if (pos->children[key[x] - 'A'] == nullptr) //If a position does not exist, return the empty vector
			{
				return totValue;
			}
			pos = pos->children[key[x] - 'A'];
		}
		return pos->values; //Return the last vector
	}
	else
	{
		for (int childCount = 0; childCount < 26; childCount++) //Iterate through all the children
		{
			//If there is a child and the string is greater than a charater...
			if (pos->children[childCount] != nullptr && key.length() > 1)
			{
				//If the character is the same as the key character, recurse with exactMatchOnly false as the character has not been switched
				if (childCount == key[0] - 'A') 
				{
					vector<ValueType> totValueHelper;
					totValueHelper = findHelper(nextKey, false, pos->children[childCount]);
					totValue.insert(totValue.end(), totValueHelper.begin(), totValueHelper.end());
				}
				//Otherwise, recurse with exactMatchOnly true as the character has been switched
				else
				{
					vector<ValueType> totValueHelper;
					totValueHelper = findHelper(nextKey, true, pos->children[childCount]);
					totValue.insert(totValue.end(), totValueHelper.begin(), totValueHelper.end());
				}
			}
			//If the string only has one character, add to the total vectors the values each legitimate child has
			else if (pos->children[childCount] != nullptr && key.length() == 1)
			{
				totValue.insert(totValue.end(), pos->children[childCount]->values.begin(), pos->children[childCount]->values.end());
			}
		}
		return totValue;
	}
}

template<typename ValueType>
vector<ValueType> Trie<ValueType>::find(const std::string& key, bool exactMatchOnly) const
{
	const Trie* pos = this;
	vector<ValueType> totValue;

	if (key == "") //Take care of empty string case
	{
		return totValue;
	}


	if (pos->children[key[0] - 'A'] == nullptr) //If the first character of the key does not have a match, return emptry vector
	{
		return totValue;
	}

	if (key.length() == 1) //Take care of one char string case
	{
		return pos->children[key[0] - 'A']->values;
	}

	pos = pos->children[key[0] - 'A'];

	if (exactMatchOnly)
	{
		return findHelper(key.substr(1), true, pos);
	}
	else
	{
		return findHelper(key.substr(1), false, pos);
	}
}
#endif // TRIE_INCLUDED
