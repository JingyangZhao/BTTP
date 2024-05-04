#include <iostream>
#include <cstdio>
#include <algorithm>
#include <vector>
#include <ctime>
#include <cmath>
#include <unordered_set>
#include <bitset>
#include <unordered_map>

using namespace std;

int n = 16;

double haversineDistance(double x1, double y1, double x2, double y2) 
{
	
	pair<double, double> coord1, coord2;
	
	coord1.first = x1, coord1.second = y1;
	coord2.first = x2, coord2.second = y2;
	
    const double R = 3959.0;  // Earth radius in miles

    double lat1 = coord1.first * M_PI / 180.0;
    double lon1 = coord1.second * M_PI / 180.0;
    double lat2 = coord2.first * M_PI / 180.0;
    double lon2 = coord2.second * M_PI / 180.0;

    double dlat = lat2 - lat1;
    double dlon = lon2 - lon1;

    double a = sin(dlat / 2) * sin(dlat / 2) + cos(lat1) * cos(lat2) * sin(dlon / 2) * sin(dlon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));

    return R * c;
}

void load(vector<vector<double> >& matrix){

    vector<vector<double> > coordinates(2*n+1, vector<double>(2+1));
    FILE *fp = NULL;
    fp = fopen("locations.txt","r");
    
	for(int i=1; i<=2*n; i++){
		fscanf(fp,"%lf %lf", &coordinates[i][1], &coordinates[i][2]);
	}

    fclose(fp);

	for(int i=1; i<=2*n; i++){
		for(int j=1; j<=2*n; j++){
			matrix[i][j] = haversineDistance(coordinates[i][1], coordinates[i][2], coordinates[j][1], coordinates[j][2]);
		}
	}

}

struct SetHash {
    size_t operator()(const std::unordered_set<int>& s) const {
        size_t hash = 0;
        for (int elem : s) {

            hash += 1 << elem;

            //hash ^= std::hash<int>{}(elem);
        }
        return hash;
    }
};

int settoint(unordered_set<int>& st){

    int ans = 0;

    for (auto& elem : st){
            ans += (1 << (elem-1));
    }
    
    return ans;
}

double cvrp(vector<vector<double> >& matrix)
{
    int s = 1; //s = 2^n;

    for(int i=1; i<=n; i++){
        s = s*2;
    }

    double ans = 0;

    unordered_map<unordered_set<int>, double, SetHash> opt; // i to opt(set(i))

    vector<unordered_set<int>> f(s+1); // i to set(i)

    vector<vector<unordered_set<int>>> S(n+1); // i to sets of size i

    for(int i=1; i<=s; i++)
    {
        int tem = 1;

        int t=0;
        for(int j=1; j<=n; j++)
        {
            if((int)(i&tem) == tem)
            {
                f[i].insert(j);
                t++;
            }

            tem = tem << 1;
        }

        S[t].push_back(f[i]);
    }

    cout<<"start:"<<endl;
    for(int i=1; i<=n; i++){
        cout<<i;

        for(auto& st: S[i]){
            if(st.size()==1){
                vector<int> lst;
                for(auto& vertex: st){
                    lst.push_back(vertex);
                }
                opt[st] = 2*matrix[0][lst[0]];
            }
            else if(st.size()==2){
                vector<int> lst;
                for(auto& vertex: st){
                    lst.push_back(vertex);
                }
                opt[st] = matrix[0][lst[0]] + matrix[lst[0]][lst[1]] + matrix[0][lst[1]];
            }
            else if(st.size()==3){
                vector<int> lst;
                for(auto& vertex: st){
                    lst.push_back(vertex);
                }
                double tem = 1e10; //, tem2 = 1e10, tem3 = 1e10;
                //sort(lst.begin(), lst.end());

                tem = min(tem, matrix[0][lst[0]] + matrix[lst[0]][lst[1]] + matrix[lst[1]][lst[2]] + matrix[0][lst[2]]);

                swap(lst[0], lst[1]);

                tem = min(tem, matrix[0][lst[0]] + matrix[lst[0]][lst[1]] + matrix[lst[1]][lst[2]] + matrix[0][lst[2]]);

                swap(lst[1], lst[2]);

                tem = min(tem, matrix[0][lst[0]] + matrix[lst[0]][lst[1]] + matrix[lst[1]][lst[2]] + matrix[0][lst[2]]);

                opt[st] = tem;
            }
            else{
                vector<int> lst;
                for(auto& vertex: st){
                    lst.push_back(vertex);
                }

                double tem1 = 1e10, tem2 = 1e10, tem3 = 1e10;
                
                for(int l1=0; l1<lst.size(); l1++){
                    st.erase(lst[l1]);

                    tem1 = min(tem1, 2*matrix[0][lst[l1]] + opt[st]);

                    for(int l2=l1+1; l2<lst.size(); l2++){
                        //if(l2 == l1) continue;

                        st.erase(lst[l2]);

                        tem2 = min(tem2, matrix[0][lst[l1]] + matrix[lst[l1]][lst[l2]] + matrix[0][lst[l2]] + opt[st]);

                        for(int l3=l2+1; l3<lst.size(); l3++){
                            //if(l3 == l1 || l3 == l2) continue;

                            st.erase(lst[l3]);

                            vector<int> lstt;
                            lstt.push_back(lst[l1]);
                            lstt.push_back(lst[l2]);
                            lstt.push_back(lst[l3]);


                            tem3 = min(tem3, matrix[0][lstt[0]] + matrix[lstt[0]][lstt[1]] + matrix[lstt[1]][lstt[2]] + matrix[0][lstt[2]] + opt[st]);

                            swap(lstt[0], lstt[1]);

                            tem3 = min(tem3, matrix[0][lstt[0]] + matrix[lstt[0]][lstt[1]] + matrix[lstt[1]][lstt[2]] + matrix[0][lstt[2]] + opt[st]);

                            swap(lstt[1], lstt[2]);

                            tem3 = min(tem3, matrix[0][lstt[0]] + matrix[lstt[0]][lstt[1]] + matrix[lstt[1]][lstt[2]] + matrix[0][lstt[2]] + opt[st]);

                            st.insert(lst[l3]);

                        }

                        st.insert(lst[l2]);
                    }

                    st.insert(lst[l1]);
                }
                double tem = min(tem1, min(tem2, tem3));
                opt[st] = tem;
                if(st.size()==n) ans = tem;
            }
        }
        cout<<": end"<<endl;
    }

    //unordered_set<int> wh;

    //for(int i=1; i<=n; i++) wh.insert(i); 

    //return opt[wh];

    return ans;
}

int main(){

    clock_t s_clock, t_clock; 

	s_clock = clock();

    vector<vector<double> > matrix(2*n+1, vector<double>(2*n+1));

    load(matrix);
	
    vector<double> ilbx(n+1), ilby(n+1);

    vector<vector<double>> newmatrix(n+1, vector<double>(n+1));

    // test for a single team x_1
    // unordered_map<int, int> mp;
    // int i=1;
    //     mp[0] = i;
    //     for(int j=1; j<=n; j++) mp[j] = n+j;

    //     for(int i=0; i<=n; i++){
    //         for(int j=0; j<=n; j++){
    //             newmatrix[i][j] = matrix[mp[i]][mp[j]];
    //         }
    //     }
    // cout<<cvrp(newmatrix);
    
    for(int i=1; i<=n; i++){

        cout<<"*********"<<i<<endl;

        unordered_map<int, int> mp;

        mp[0] = i;
        for(int j=1; j<=n; j++) mp[j] = n+j;

        for(int j=0; j<=n; j++){
            for(int k=0; k<=n; k++){
                newmatrix[j][k] = matrix[mp[j]][mp[k]];
            }
        }

        ilbx[i] = cvrp(newmatrix);

        mp[0] = i+n;
        for(int j=1; j<=n; j++) mp[j] = j;

        for(int j=0; j<=n; j++){
            for(int k=0; k<=n; k++){
                newmatrix[j][k] = matrix[mp[j]][mp[k]];
            }
        }

        ilby[i] = cvrp(newmatrix);
    }

    t_clock = clock();

    double ans = 0.0;

    for(int i=1; i<=n; i++){
        ans += ilbx[i];
        ans += ilby[i];
    }

    printf("%.4lf\n", ans);

    printf("The running time is:%lf\n", double(t_clock-s_clock)/CLOCKS_PER_SEC);

    system("pause");

    return 0;
}