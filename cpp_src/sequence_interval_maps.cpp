#include <iostream>
#include <boost/icl/interval_map.hpp>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;
using namespace boost;
using namespace boost::icl;

void boost_intervals(vector<int> start, vector<int> end, vector<float> p_vector)
{
    interval_map<int, float> sequence_interval_map;
    
    int tmp_start, tmp_end, i;
    float tmp_psi;
    for(i=0; i<(int) start.size(); ++i){
        // get all values that are together in a row
        tmp_start = start.at(i);
        tmp_end = end.at(i);
        tmp_psi = p_vector.at(i);
        
        sequence_interval_map += // element addition can also be done via operator +=
          make_pair( 
            interval<int>::right_open(
              tmp_start, 
              tmp_end), 
            tmp_psi);
        // cout<<(*i)<<std::endl; 
    }


    ofstream output_file("output.txt");
    if(output_file.is_open())
    {
        interval_map<int, float>::iterator it = sequence_interval_map.begin();
        while(it != sequence_interval_map.end())
        {
            interval<int>::type sequence_interval = it->first;
            float psi_value = (*it++).second;
            output_file << first(sequence_interval) << "\t" << upper(sequence_interval) << "\t" << psi_value << endl;
        }
        output_file.close();
    }else
    {
        cout << "file did not open\n";
    }
}


int main(int argc, char *argv[])
{
    // define variables
    string line;
    int counter;
    ifstream csvFile ("input.txt");
    vector<int> starts;
    vector<int> ends;
    vector<float> psi_vector;

    // check if it opened
    if (csvFile.is_open())
    {
        // loop until no more lines
        while (csvFile.good())
        {
            getline(csvFile, line);
            tokenizer<escaped_list_separator<char> > tok(line);
            counter = 0;  // set the token counter back to zero
            for(tokenizer<escaped_list_separator<char> >::iterator beg=tok.begin(); beg!=tok.end();++beg){
                if(counter==0)
                {
                    // start interval column
                    starts.push_back(atoi((*beg).c_str()));
                }else if(counter==1)
                {
                    // end interval column
                    ends.push_back(atoi((*beg).c_str()));
                }else if(counter==2)
                {
                    // psi column
                    psi_vector.push_back(atof((*beg).c_str()));
                }
                counter++;
            }
        }
        
        // test print code
        // vector<int>::const_iterator i;
        // for(i=starts.begin(); i!=starts.end(); ++i){
        //    cout<<(*i)<<std::endl; 
        //}

    }
        
    boost_intervals(starts, ends, psi_vector);
    return 0;
}
