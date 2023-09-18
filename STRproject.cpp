#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
using namespace std;

//Apply the algorithm to tell how similar how two strings are by calcualting the minimum changes required to convert one string to another
int levenshteinDistance(const string& str1, const string& str2)
{
    int first = str1.length();
    int second = str2.length();

    vector<vector<int>> d(first + 1, vector<int>(second + 1)); //2D vector that stores distance between substrings

    for (int i = 0; i <= first; ++i)
        d[i][0] = i;
    //Initializing
    for (int j = 0; j <= second; ++j)
        d[0][j] = j;

    for (int i = 1; i <= first; ++i) {
        for (int j = 1; j <= second; ++j) {
            if (str1[i - 1] == str2[j - 1])
                d[i][j] = d[i - 1][j - 1];
            else
                d[i][j] = 1 + min({ d[i - 1][j], d[i][j - 1], d[i - 1][j - 1] });
        }
    }

    return d[first][second];
}

//Gets the index of the starting nucleotide that matches with its reference
int getPosOfSimilariity(const string& target, const string& test)
{
    int targetLen = target.length();
    int testLen = test.length();
    int pos = -1;
    for (int i = 0; i <= targetLen - testLen; i++)
    {
        string substr = target.substr(i, testLen); //Cutting the sequence so the lengths are the same
        int distance = levenshteinDistance(substr, test);
        if (distance <= 1)
        {
            pos = i;
            break;
        }
    }
    return pos;
}

//Inputting the starting index and the size, it returns the range of the fragmented sequence from the FASTQ to match against the reference
vector<pair<int, int>> getIndicesRange(int index, int length)
{
    vector<pair<int, int>> lengthRange;
    lengthRange.push_back(make_pair(index, index + length - 1));
    return lengthRange;
}

//Returns the position of where a mismatch/SNP happen
int mismatchPos(string target, string test)
{
    for (int i = 0; i < target.size(); i++)
    {
        if (target.at(i) != test.at(i))
        {
            return static_cast<int>(i);
        }
    }
}

//Returns the edit distance number to a corresponding cut sequece to the reference
int hasExactMatch(const string& target, const string& test)
{
    int targetLen = target.length();
    int testLen = test.length();


    for (int i = 0; i <= targetLen - testLen; i++)
    {
        string substr = target.substr(i, testLen);
        int distance = levenshteinDistance(substr, test);
        if (distance == 0) //When exactly the same
            return distance;
    }
    for (int i = 0; i <= targetLen - testLen; i++)
    {
        string substr = target.substr(i, testLen);
        int distance = levenshteinDistance(substr, test);
        if (distance == 1) //When there is one off
            return distance;
    }
    return -1;
}

//Final step for displaying purposes
//Adjusts the name of reverse complementary to its original and adds the count of the number of times 
//either the original or reverse complementary of the allele has matched to its reference
void printCounter(vector<pair<string, int>>& counter, string target)
{
    int index = target.find("_R");
    if (index != string::npos)
    {
        target = target.substr(0, index); //cuts off the _R so the variable names stay the same -> counts for original and reverse complement are combined
    }

    bool lastpos = true;
    if (counter.size() == 0)
    {
        counter.push_back(make_pair(target, 1));
    }
    else
    {
        for (int i = 0; i < counter.size(); i++)
        {
            if (counter[i].first == target)
            {
                counter[i].second += 1;
                lastpos = false;
            }
        }
        if (lastpos == true)
        {
            counter.push_back(make_pair(target, 1));
        }
    }

}

//The main function that determines whether the allele from FASTQ closely resembles with reference
bool matchSimilarity(string target, string allele, string name, int allelesize, string actualname)
//matchSimilarity(DNA[i], allele string ex: v.second[i], label[i], v.second[i].size(), v.first[i])
{
    string fastqSeq = target;
    string labelName = name;
    int  firstPos = -1;
    int  secondPos = -1;
    string  cutTarget = "";
    int exactDis = hasExactMatch(target, allele);
    //Case when there is exactly a match
    if (exactDis == 0)
    {
        cout << "Label: " << name << endl;
        cout << "The allele is " << actualname << endl;
        return true;
    }
    //Case when the edit distance is 1
    if (exactDis == 1)
    {
        int posNum = getPosOfSimilariity(target, allele); //Gets the index of the starting nucleotide
        vector<pair<int, int>> SizeSNP = getIndicesRange(posNum, allele.size());
        for (const auto& indices : SizeSNP)
        {
            firstPos = indices.first;
            secondPos = indices.second; //gets the vector pairs
        }
        cutTarget = target.substr(firstPos, secondPos - firstPos + 1); //Gets the sequence frame of FASTQ that matches with the reference in length 
        int mismatchp = mismatchPos(cutTarget, allele); //The position of where the SNP happened 
        if ((mismatchp > 8) && (mismatchp < (cutTarget.size() - 10))) //Mismatch has to occur between first 9 to last 9 nucleotides to count
        {
            cout << "Label: " << name << endl;
            cout << "The allele is " << actualname << " (with one substitution)" << endl;
            cout << "First line is from fastq, second line is from reference allele " << endl;
            cout << cutTarget << endl;
            cout << allele << endl;
            return true;
        }
    }
    return false;

}



int main() {

    //Loading in the possible alleles and its reverse complement as strings into a vector pair of name and its nucleotides
    //The allele database has been sorted from longest to shortest
    vector <pair<string, string>> v;

    string AMELX = "GGAAGCTGGTGGTAGGAACTGTAAAATCAGGACCACTTGAGAAACATCTGGGATAAAGAATCAACACACTATTCTTTACAGAGCCCAGGGCATTGTTAACGCAAACAATGGTCAAAATTA";
    v.push_back(make_pair("AMELX", AMELX));
    string AMELX_R = "TAATTTTGACCATTGTTTGCGTTAACAATGCCCTGGGCTCTGTAAAGAATAGTGTGTTGATTCTTTATCCCAGATGTTTCTCAAGTGGTCCTGATTTTACAGTTCCTACCACCAGCTTCC";
    v.push_back(make_pair("AMELX_R", AMELX_R));
    string AMELY = "TAGGAACTGTAAAATTGGGACCACTTGAGAAACCACTTTATTTGGGATGAAGAATCCACCCACTATTCTTTACAGAGCCCAGGGGACTGCTAATGCAAACAGTGATCAAAATTAGTAAAG";
    v.push_back(make_pair("AMELY", AMELY));
    string AMELY_R = "CTTTACTAATTTTGATCACTGTTTGCATTAGCAGTCCCCTGGGCTCTGTAAAGAATAGTGGGTGGATTCTTCATCCCAAATAAAGTGGTTTCTCAAGTGGTCCCAATTTTACAGTTCCTA";
    v.push_back(make_pair("AMELY_R", AMELY_R));
    string D3S1744_21 = "ACCTGTCTATCTATCTATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCATCTATCTATAG";
    v.push_back(make_pair("D3S1744_21", D3S1744_21));
    string D3S1744_21_R = "CTATAGATAGATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATATAGATAGATAGACAGGT";
    v.push_back(make_pair("D3S1744_21_R", D3S1744_21_R));
    string D3S1744_20 = "ACCTGTCTATCTATCTATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCACTTATCTATAG";
    v.push_back(make_pair("D3S1744_20", D3S1744_20));
    string D3S1744_20_R = "CTATAGATAAGTGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATATAGATAGATAGACAGGT";
    v.push_back(make_pair("D3S1744_20_R", D3S1744_20_R));
    string D3S1744_19 = "ACCTGTCTATCTATCTATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCACTTATCTATAG";
    v.push_back(make_pair("D3S1744_19", D3S1744_19));
    string D3S1744_19_R = "CTATAGATAAGTGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATATAGATAGATAGACAGGT";
    v.push_back(make_pair("D3S1744_19_R", D3S1744_19_R));
    string D3S1744_18 = "ACCTGTCTATCTATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCATCTATCTATAG";
    v.push_back(make_pair("D3S1744_18", D3S1744_18));
    string D3S1744_18_R = "CTATAGATAGATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATATAGATAGACAGGT";
    v.push_back(make_pair("D3S1744_18_R", D3S1744_18_R));
    string D3S1744_17 = "ACCTGTCTATCTATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCATCTATCTATCTATAG";
    v.push_back(make_pair("D3S1744_17", D3S1744_17));
    string D3S1744_17_R = "CTATAGATAGATAGATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATATAGATAGACAGGT";
    v.push_back(make_pair("D3S1744_17_R", D3S1744_17_R));
    string D21S1437_17 = "AGGAGGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGACA";
    v.push_back(make_pair("D21S1437_17", D21S1437_17));
    string D21S1437_17_R = "TGTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCCTCCT";
    v.push_back(make_pair("D21S1437_17_R", D21S1437_17_R));
    string D12S1090_16_3 = "ATATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATGATAGATAGATAGAGCG";
    v.push_back(make_pair("D12S1090_16_3", D12S1090_16_3));
    string D12S1090_16_3_R = "CGCTCTATCTATCTATCATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATAT";
    v.push_back(make_pair("D12S1090_16_3_R", D12S1090_16_3_R));
    string D3S1744_16 = "ACCTGTCTATCTATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCATCTATCTATAG";
    v.push_back(make_pair("D3S1744_16", D3S1744_16));
    string D3S1744_16_R = "CTATAGATAGATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATATAGATAGACAGGT";
    v.push_back(make_pair("D3S1744_16_R", D3S1744_16_R));
    string D21S1437_16 = "AGGAGGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGACA";
    v.push_back(make_pair("D21S1437_16", D21S1437_16));
    string D21S1437_16_R = "TGTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCCTCCT";
    v.push_back(make_pair("D21S1437_16_R", D21S1437_16_R));
    string D12S1090_15_3 = "ATATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATGATAGATAGATAGAGCG";
    v.push_back(make_pair("D12S1090_15_3", D12S1090_15_3));
    string D12S1090_15_3_R = "CGCTCTATCTATCTATCATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATAT";
    v.push_back(make_pair("D12S1090_15_3_R", D12S1090_15_3_R));
    string D3S1744_15 = "ACCTGTCTATCTATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCATCTATCTATAG";
    v.push_back(make_pair("D3S1744_15", D3S1744_15));
    string D3S1744_15_R = "CTATAGATAGATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATATAGATAGACAGGT";
    v.push_back(make_pair("D3S1744_15_R", D3S1744_15_R));
    string D21S1437_15 = "AGGAGGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGACA";
    v.push_back(make_pair("D21S1437_15", D21S1437_15));
    string D21S1437_15_R = "TGTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCCTCCT";
    v.push_back(make_pair("D21S1437_15_R", D21S1437_15_R));
    string D12S1090_15 = "ATATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAGCG";
    v.push_back(make_pair("D12S1090_15", D12S1090_15));
    string D12S1090_15_R = "CGCTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATAT";
    v.push_back(make_pair("D12S1090_15_R", D12S1090_15_R));
    string D4S2366_15 = "AATGGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATTGATTGATAGACGAT";
    v.push_back(make_pair("D4S2366_15", D4S2366_15));
    string D4S2366_15_R = "ATCGTCTATCAATCAATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCATT";
    v.push_back(make_pair("D4S2366_15_R", D4S2366_15_R));
    string D12S1090_14_3 = "ATATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATGATAGATAGATAGAGCG";
    v.push_back(make_pair("D12S1090_14_3", D12S1090_14_3));
    string D12S1090_14_3_R = "CGCTCTATCTATCTATCATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATAT";
    v.push_back(make_pair("D12S1090_14_3_R", D12S1090_14_3_R));
    string D3S1744_14 = "ACCTGTCTATCTATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCATCTATCTATAG";
    v.push_back(make_pair("D3S1744_14", D3S1744_14));
    string D3S1744_14_R = "CTATAGATAGATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATATAGATAGACAGGT";
    v.push_back(make_pair("D3S1744_14_R", D3S1744_14_R));
    string D21S1437_14 = "AGGAGGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGACA";
    v.push_back(make_pair("D21S1437_14", D21S1437_14));
    string D21S1437_14_R = "TGTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCCTCCT";
    v.push_back(make_pair("D21S1437_14_R", D21S1437_14_R));
    string D12S1090_14 = "ATATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAGCG";
    v.push_back(make_pair("D12S1090_14", D12S1090_14));
    string D12S1090_14_R = "CGCTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATAT";
    v.push_back(make_pair("D12S1090_14_R", D12S1090_14_R));
    string D14S608_14 = "TATACTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTGT";
    v.push_back(make_pair("D14S608_14", D14S608_14));
    string D14S608_14_R = "ACAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAGTATA";
    v.push_back(make_pair("D14S608_14_R", D14S608_14_R));
    string D4S2366_14 = "AATGGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATTGATTGATAGACGAT";
    v.push_back(make_pair("D4S2366_14", D4S2366_14));
    string D4S2366_14_R = "ATCGTCTATCAATCAATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCATT";
    v.push_back(make_pair("D4S2366_14_R", D4S2366_14_R));
    string D18S536_14 = "ATAGGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGA";
    v.push_back(make_pair("D18S536_14", D18S536_14));
    string D18S536_14_R = "TCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCTAT";
    v.push_back(make_pair("D18S536_14_R", D18S536_14_R));
    string D12S1090_13_3 = "ATATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATGATAGATAGATAGAGCG";
    v.push_back(make_pair("D12S1090_13_3", D12S1090_13_3));
    string D12S1090_13_3_R = "CGCTCTATCTATCTATCATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATAT";
    v.push_back(make_pair("D12S1090_13_3_R", D12S1090_13_3_R));
    string D3S1744_13 = "ACCTGTCTATCTATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCATCTATCTATAG";
    v.push_back(make_pair("D3S1744_13", D3S1744_13));
    string D3S1744_13_R = "CTATAGATAGATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATATAGATAGACAGGT";
    v.push_back(make_pair("D3S1744_13_R", D3S1744_13_R));
    string D21S1437_13 = "AGGAGGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGACA";
    v.push_back(make_pair("D21S1437_13", D21S1437_13));
    string D21S1437_13_R = "TGTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCCTCCT";
    v.push_back(make_pair("D21S1437_13_R", D21S1437_13_R));
    string D12S1090_13 = "ATATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAGCG";
    v.push_back(make_pair("D12S1090_13", D12S1090_13));
    string D12S1090_13_R = "CGCTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATAT";
    v.push_back(make_pair("D12S1090_13_R", D12S1090_13_R));
    string D14S608_13 = "TATACTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTGT";
    v.push_back(make_pair("D14S608_13", D14S608_13));
    string D14S608_13_R = "ACAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAGTATA";
    v.push_back(make_pair("D14S608_13_R", D14S608_13_R));
    string D4S2366_13 = "AATGGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATTGATTGATAGACGAT";
    v.push_back(make_pair("D4S2366_13", D4S2366_13));
    string D4S2366_13_R = "ATCGTCTATCAATCAATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCATT";
    v.push_back(make_pair("D4S2366_13_R", D4S2366_13_R));
    string D18S536_13 = "ATAGGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGA";
    v.push_back(make_pair("D18S536_13", D18S536_13));
    string D18S536_13_R = "TCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCTAT";
    v.push_back(make_pair("D18S536_13_R", D18S536_13_R));
    string D13S765_13 = "ATCATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATGGAAT";
    v.push_back(make_pair("D13S765_13", D13S765_13));
    string D13S765_13_R = "ATTCCATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATGAT";
    v.push_back(make_pair("D13S765_13_R", D13S765_13_R));
    string D12S1090_12_3 = "ATATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATGATAGATAGATAGAGCG";
    v.push_back(make_pair("D12S1090_12_3", D12S1090_12_3));
    string D12S1090_12_3_R = "CGCTCTATCTATCTATCATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATAT";
    v.push_back(make_pair("D12S1090_12_3_R", D12S1090_12_3_R));
    string D21S1437_12 = "AGGAGGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGACA";
    v.push_back(make_pair("D21S1437_12", D21S1437_12));
    string D21S1437_12_R = "TGTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCCTCCT";
    v.push_back(make_pair("D21S1437_12_R", D21S1437_12_R));
    string D12S1090_12 = "ATATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAGCG";
    v.push_back(make_pair("D12S1090_12", D12S1090_12));
    string D12S1090_12_R = "CGCTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATAT";
    v.push_back(make_pair("D12S1090_12_R", D12S1090_12_R));
    string D14S608_12 = "TATACTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTGT";
    v.push_back(make_pair("D14S608_12", D14S608_12));
    string D14S608_12_R = "ACAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAGTATA";
    v.push_back(make_pair("D14S608_12_R", D14S608_12_R));
    string D4S2366_12 = "AATGGATAGATAGATAGATAGATAGATAGATAGATAGATAGATTGATTGATAGACGAT";
    v.push_back(make_pair("D4S2366_12", D4S2366_12));
    string D4S2366_12_R = "ATCGTCTATCAATCAATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCATT";
    v.push_back(make_pair("D4S2366_12_R", D4S2366_12_R));
    string D18S536_12 = "ATAGGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGA";
    v.push_back(make_pair("D18S536_12", D18S536_12));
    string D18S536_12_R = "TCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCTAT";
    v.push_back(make_pair("D18S536_12_R", D18S536_12_R));
    string D13S765_12 = "ATCATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATGGAAT";
    v.push_back(make_pair("D13S765_12", D13S765_12));
    string D13S765_12_R = "ATTCCATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATGAT";
    v.push_back(make_pair("D13S765_12_R", D13S765_12_R));
    string D21S1437_11 = "AGGAGGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGACA";
    v.push_back(make_pair("D21S1437_11", D21S1437_11));
    string D21S1437_11_R = "TGTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCCTCCT";
    v.push_back(make_pair("D21S1437_11_R", D21S1437_11_R));
    string D12S1090_11 = "ATATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAGCG";
    v.push_back(make_pair("D12S1090_11", D12S1090_11));
    string D12S1090_11_R = "CGCTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATAT";
    v.push_back(make_pair("D12S1090_11_R", D12S1090_11_R));
    string D14S608_11 = "TATACTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTGT";
    v.push_back(make_pair("D14S608_11", D14S608_11));
    string D14S608_11_R = "ACAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAGTATA";
    v.push_back(make_pair("D14S608_11_R", D14S608_11_R));
    string D4S2366_11 = "AATGGATAGATAGATAGATAGATAGATAGATAGATAGATTGATTGATAGACGAT";
    v.push_back(make_pair("D4S2366_11", D4S2366_11));
    string D4S2366_11_R = "ATCGTCTATCAATCAATCTATCTATCTATCTATCTATCTATCTATCTATCCATT";
    v.push_back(make_pair("D4S2366_11_R", D4S2366_11_R));
    string D18S536_11 = "ATAGGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGA";
    v.push_back(make_pair("D18S536_11", D18S536_11));
    string D18S536_11_R = "TCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCTAT";
    v.push_back(make_pair("D18S536_11_R", D18S536_11_R));
    string D13S765_11 = "ATCATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATGGAAT";
    v.push_back(make_pair("D13S765_11", D13S765_11));
    string D13S765_11_R = "ATTCCATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATGAT";
    v.push_back(make_pair("D13S765_11_R", D13S765_11_R));
    string D12S1090_10_3 = "ATATAGATAGATAGATAGATAGATAGATAGATGATAGATAGATAGATAGAGCG";
    v.push_back(make_pair("D12S1090_10_3", D12S1090_10_3));
    string D12S1090_10_3_R = "CGCTCTATCTATCTATCTATCATCTATCTATCTATCTATCTATCTATCTATAT";
    v.push_back(make_pair("D12S1090_10_3_R", D12S1090_10_3_R));
    string D21S1437_10 = "AGGAGGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGACA";
    v.push_back(make_pair("D21S1437_10", D21S1437_10));
    string D21S1437_10_R = "TGTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCCTCCT";
    v.push_back(make_pair("D21S1437_10_R", D21S1437_10_R));
    string D12S1090_10 = "ATATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAGCG";
    v.push_back(make_pair("D12S1090_10", D12S1090_10));
    string D12S1090_10_R = "CGCTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATAT";
    v.push_back(make_pair("D12S1090_10_R", D12S1090_10_R));
    string D14S608_10 = "TATACTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTGT";
    v.push_back(make_pair("D14S608_10", D14S608_10));
    string D14S608_10_R = "ACAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAGTATA";
    v.push_back(make_pair("D14S608_10_R", D14S608_10_R));
    string D4S2366_10 = "AATGGATAGATAGATAGATAGATAGATAGATAGATTGATTGATAGACGAT";
    v.push_back(make_pair("D4S2366_10", D4S2366_10));
    string D4S2366_10_R = "ATCGTCTATCAATCAATCTATCTATCTATCTATCTATCTATCTATCCATT";
    v.push_back(make_pair("D4S2366_10_R", D4S2366_10_R));
    string D18S536_10 = "ATAGGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGA";
    v.push_back(make_pair("D18S536_10", D18S536_10));
    string D18S536_10_R = "TCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCTAT";
    v.push_back(make_pair("D18S536_10_R", D18S536_10_R));
    string D13S765_10 = "ATCATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATGGAAT";
    v.push_back(make_pair("D13S765_10", D13S765_10));
    string D13S765_10_R = "ATTCCATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATGAT";
    v.push_back(make_pair("D13S765_10_R", D13S765_10_R));
    string D21S1437_9 = "AGGAGGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGACA";
    v.push_back(make_pair("D21S1437_9", D21S1437_9));
    string D21S1437_9_R = "TGTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCCTCCT";
    v.push_back(make_pair("D21S1437_9_R", D21S1437_9_R));
    string D12S1090_9 = "ATATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAGCG";
    v.push_back(make_pair("D12S1090_9", D12S1090_9));
    string D12S1090_9_R = "CGCTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATAT";
    v.push_back(make_pair("D12S1090_9_R", D12S1090_9_R));
    string D14S608_9 = "TATACTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTGT";
    v.push_back(make_pair("D14S608_9", D14S608_9));
    string D14S608_9_R = "ACAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAGTATA";
    v.push_back(make_pair("D14S608_9_R", D14S608_9_R));
    string D4S2366_9 = "AATGGATAGATAGATAGATAGATAGATAGATTGATTGATAGACGAT";
    v.push_back(make_pair("D4S2366_9", D4S2366_9));
    string D4S2366_9_R = "ATCGTCTATCAATCAATCTATCTATCTATCTATCTATCTATCCATT";
    v.push_back(make_pair("D4S2366_9_R", D4S2366_9_R));
    string D18S536_9 = "ATAGGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGA";
    v.push_back(make_pair("D18S536_9", D18S536_9));
    string D18S536_9_R = "TCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCCTAT";
    v.push_back(make_pair("D18S536_9_R", D18S536_9_R));
    string D13S765_9 = "ATCATAGATAGATAGATAGATAGATAGATAGATAGATAGATGGAAT";
    v.push_back(make_pair("D13S765_9", D13S765_9));
    string D13S765_9_R = "ATTCCATCTATCTATCTATCTATCTATCTATCTATCTATCTATGAT";
    v.push_back(make_pair("D13S765_9_R", D13S765_9_R));
    string D14S608_8 = "TATACTCTATCTATCTATCTATCTATCTATCTATCTATCTGT";
    v.push_back(make_pair("D14S608_8", D14S608_8));
    string D14S608_8_R = "ACAGATAGATAGATAGATAGATAGATAGATAGATAGAGTATA";
    v.push_back(make_pair("D14S608_8_R", D14S608_8_R));
    string D18S536_8 = "ATAGGATAGATAGATAGATAGATAGATAGATAGATAGACAGA";
    v.push_back(make_pair("D18S536_8", D18S536_8));
    string D18S536_8_R = "TCTGTCTATCTATCTATCTATCTATCTATCTATCTATCCTAT";
    v.push_back(make_pair("D18S536_8_R", D18S536_8_R));
    string D13S765_8 = "ATCATAGATAGATAGATAGATAGATAGATAGATAGATGGAAT";
    v.push_back(make_pair("D13S765_8", D13S765_8));
    string D13S765_8_R = "ATTCCATCTATCTATCTATCTATCTATCTATCTATCTATGAT";
    v.push_back(make_pair("D13S765_8_R", D13S765_8_R));
    string D21S1437_7 = "AGGAGGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGACA";
    v.push_back(make_pair("D21S1437_7", D21S1437_7));
    string D21S1437_7_R = "TGTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCCTCCT";
    v.push_back(make_pair("D21S1437_7_R", D21S1437_7_R));
    string D12S1090_7 = "ATATAGATAGATAGATAGATAGATAGATAGATAGAGCG";
    v.push_back(make_pair("D12S1090_7", D12S1090_7));
    string D12S1090_7_R = "CGCTCTATCTATCTATCTATCTATCTATCTATCTATAT";
    v.push_back(make_pair("D12S1090_7_R", D12S1090_7_R));
    string D14S608_7 = "TATACTCTATCTATCTATCTATCTATCTATCTATCTGT";
    v.push_back(make_pair("D14S608_7", D14S608_7));
    string D14S608_7_R = "ACAGATAGATAGATAGATAGATAGATAGATAGAGTATA";
    v.push_back(make_pair("D14S608_7_R", D14S608_7_R));
    string D13S765_7 = "ATCATAGATAGATAGATAGATAGATAGATAGATGGAAT";
    v.push_back(make_pair("D13S765_7", D13S765_7));
    string D13S765_7_R = "ATTCCATCTATCTATCTATCTATCTATCTATCTATGAT";
    v.push_back(make_pair("D13S765_7_R", D13S765_7_R));
    string D12S1090_6_R = "CGCTCTATCTATCTATCTATCTATCTATCTATAT";
    v.push_back(make_pair("D12S1090_6_R", D12S1090_6_R));
    string D12S1090_6 = "ATATAGATAGATAGATAGATAGATAGATAGAGCG";
    v.push_back(make_pair("D12S1090_6", D12S1090_6));
    string D14S608_6 = "TATACTCTATCTATCTATCTATCTATCTATCTGT";
    v.push_back(make_pair("D14S608_6", D14S608_6));
    string D14S608_6_R = "ACAGATAGATAGATAGATAGATAGATAGAGTATA";
    v.push_back(make_pair("D14S608_6_R", D14S608_6_R));
    string D14S608_5 = "TATACTCTATCTATCTATCTATCTATCTGT";
    v.push_back(make_pair("D14S608_5", D14S608_5));
    string D14S608_5_R = "ACAGATAGATAGATAGATAGATAGAGTATA";
    v.push_back(make_pair("D14S608_5_R", D14S608_5_R));


    ////////////////////////////////////////////////////////////


    //Importing FASTQ files and inserting only the labels and sequences in separate vectors
    string input = "";
    cout << "Please enter the EXACT file name:";
    cin >> input;

    ifstream inputFile(input);
    string line;
    vector<string> DNA;
    vector<string> label;
    vector<pair<string, int>> c; //Vector pair that tracks the total counts of each allele when the program finishes running
    bool isnext = false; //This check ensures that only the line followed by the label line gets inputted into the sequence vector
    if (!inputFile.is_open())
    {
        cerr << "Unable to open or locate the file." << endl;
    }

    int counter = 0; //The counter check ensures that only the label line is inserted into the vector
    while (getline(inputFile, line)) //Looping through each line in the FASTQ and inputting only the label line and sequence line in their respective vectors
    {
        if (line[0] == '@') //When a label line or quality line that starts with ASCII @ get detected
        {
            if (counter == 0)
            {
                counter++;
                label.push_back(line);
                isnext = true;
            }
            else if (counter == 1)
                //The code only goes here when the current line is the label line and the previous inserted line is the quality line
            {
                label.pop_back(); //Pop the quality line from the vector
                counter = 0; //Reset the counter variable to reset the process
                label.push_back(line);//Insert the label line to the vector
                isnext = true;
            }
        }
        else if (isnext == true)
        {
            counter = 0;
            if ((line[0] == 'G' || line[0] == 'C' || line[0] == 'A' || line[0] == 'T')
                && (line[1] == 'G' || line[1] == 'C' || line[1] == 'A' || line[1] == 'T')
                && (line[2] == 'G' || line[2] == 'C' || line[2] == 'A' || line[2] == 'T')) //When there is a sequencing line
            {
                DNA.push_back(line);
                isnext = false;
            }
        }

    }

    cout << "==========================================================" << endl;

    for (int i = 0; i < label.size(); i++) //Looping through each label/sequence 
    {
        bool ischecked = false; //Once a line is checked, there is no need to check again
        string labelName = label[i];
        string seqTarget = DNA[i];
        int targetSize = seqTarget.length(); //Gets the DNA sequence of a label
        for (int i = 0; i < v.size(); i++)
        {
            if (ischecked == false)
            {
                if (matchSimilarity(seqTarget, v[i].second, labelName, v[i].second.size(), v[i].first) == true) //If similar
                {
                    printCounter(c, v[i].first);//Adjust any _R strings and add to paired vector for results
                    ischecked = true;
                }
            }
        }
    }

    cout << " " << endl;
    cout << "======================================================================" << endl;
    cout << " " << endl;

    if (c.size() == 0)
    {
        cout << "Nothing!" << endl;
    }
    else
    {
        for (int i = 0; i < c.size(); i++)
        {
            cout << "There are " << c[i].second << " " << c[i].first << endl;
        }
    }

    inputFile.close();
    cout << "Program finished running" << endl;

}