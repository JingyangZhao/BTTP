#include <iostream>
#include <cstdio>
#include <algorithm>
#include <vector>
#include <ctime>
#include <cmath>

using namespace std;

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
//for showing the shcedule
void show_table(int n, vector<vector<int> > &day)
{
	printf("/**************table****************/\n");
	for(int i=1; i<=2*n; i++,printf("\n"))
		for(int j=1; j<=2*n; j++)	
			printf("%d ", day[i][j]);
}

//check for the feasibility of our schedule
int check(int n, vector<vector<int> > &day)
{
	//check for whether there is a missed game  
	vector<vector<int> > c(2*n+1, vector<int>(2*n+1));
	for(int i=1; i<=n; i++) for(int j=n+1; j<=2*n; j++) c[i][j]=200;
	for(int i=n+1; i<=2*n; i++) for(int j=1; j<=n; j++) c[i][j]=200;
	
	for(int i=1; i<=n; i++)
	{
		for(int j=1; j<=2*n; j++)
		{
			if(day[i][j]>0) c[i][day[i][j]] += (-100+1);
			else c[i][-day[i][j]] += (-100-1);
		}
	}
	
	for(int i=n+1; i<=2*n; i++)
	{
		for(int j=1; j<=2*n; j++)
		{
			if(day[i][j]>0) c[i][day[i][j]] += (-100+1);
			else c[i][-day[i][j]] += (-100-1);
		}
	}
	
	for(int i=1; i<=n; i++) for(int j=1+n; j<=n+n; j++) if(c[i][j]!=0){ printf("missing: %d %d\n", i, j); return -1; }
	for(int i=n+1; i<=2*n; i++) for(int j=1; j<=n; j++) if(c[i][j]!=0){ printf("missing: %d %d\n", i, j); return -1; }
	
	//check for the condition of BTTP-3: bounded-by-3
	for(int i=1; i<=2*n; i++)
	{
		for(int d=1; d<=2*n-3; d++)
		{
			if(day[i][d]>0 && day[i][d+1]>0 && day[i][d+2]>0 && day[i][d+3]>0)
			{
				printf("bounded by-3: team:%d day:%d\n", i, d);
				return -1;
			}
			if(day[i][d]<0 && day[i][d+1]<0 && day[i][d+2]<0 && day[i][d+3]<0)
			{
				printf("bounded by-3: team:%d day:%d\n", i, d);
				return -1;
			}
		}
	}
	
	//check for no-repeat
	for(int i=1; i<=2*n; i++) 
		for (int d=1; d<=2*n-1; d++)
			if(day[i][d]==day[i][d+1] || day[i][d]== -day[i][d+1])
			{
				printf("repeat: team:%d day:%d\n", i, d);
				return -1;			
			}
			
	return 0;
} 

//calculate the total distance travelled by all teams
//according the schedule: &day
double dist(int n, vector<vector<int> > &day, vector<vector<double> > &D)
{
	double res=0;
	for(int i=1; i<=2*n; i++)
	{
		int flag=0;
		for(int d=1; d<=2*n; d++)
			if(d==2*n && day[i][d]>0)
				if(flag==0) res += 2*D[i][day[i][d]];	
				else res += D[day[i][d-1]][day[i][d]]+D[day[i][d]][i];
			else if(day[i][d]<0)
				if(flag==0) continue; 
				else res += D[i][day[i][d-1]], flag=0;
			else
				if(flag==0) res += D[i][day[i][d]], flag=1;
				else res += D[day[i][d-1]][day[i][d]];		
	}
	
	return res;
}
void normalSupergame(int above, int below, vector<int> &x, vector<int> &y, vector<vector<int> > &table, int d, int f)
{
	int n=3;
	for(int i=1; i<=2*n; i++)
	{
		if(i%6<=3 && i%6>=1) 
		{
			int opponent_1 = y[ i%6 + 3*((i-1)/6) + 3*below-3 ]; 
			table[x[1 + 3*above-3]][i + d-1] = f*opponent_1;
			table[opponent_1][i + d-1] = -f*x[1 + 3*above-3];
			
			for(int j=2; j<=n; j++)
			{
				int d_j = i + (j-1)+3*( (j-1)/3 );
				if(d_j>2*n) d_j -= 2*n;
				
				table[x[j + 3*above-3]][d_j + d-1] = f*opponent_1;
				table[opponent_1][d_j + d-1] = -f*x[j + 3*above-3];
			}
		}
		else 
		{
			int tem = (i%6==0)? 3:(i%6-3);
			int opponent_1 = y[  tem + 3*((i-1)/6) + 3*below-3];
			table[x[1 + 3*above-3]][i + d-1] = -f*opponent_1; // play home game
			table[opponent_1][i + d-1] = f*x[1 + 3*above-3];	
			for(int j=2; j<=n; j++)
			{
				int d_j = i + (j-1)+3*( (j-1)/3 );
				if(d_j>2*n) d_j -= 2*n;
				
				table[x[j + 3*above-3]][d_j + d-1] = -f*opponent_1;
				table[opponent_1][d_j + d-1] = f*x[j + 3*above-3];
			}
		}
	}
}
void leftSupergame(int n, int above, int below, vector<int> &x, vector<int> &y, vector<vector<int> > &table, int d, int f)
{
	//x[3*above-2], x[3*above-2], x[3*above-1], x[n] ------- y[3*below-2], y[3*below-2], y[3*below-1], y[n]
	int X[4]={x[3*above-2], x[3*above-1], x[3*above], x[n]};
	int Y[4]={y[3*below-2], y[3*below-1], y[3*below], y[n]};
	
	if(1){
		int Delta[3] = {1, 2, 3};
		int F[3] = {f, f, -f};
		for(int i=1; i<=3; i++){
			int delta = Delta[i-1];
			for(int j=0; j<4; j++){
				table[X[j]][d] = F[i-1]*Y[(j+delta)%4];
				table[Y[(j+delta)%4]][d] = -F[i-1]*X[j];
			}
			d++;
		}
	}
	if(1){
		int Delta[3] = {2, 3, 1};
		int F[3] = {-f, f, -f};
		for(int i=1; i<=3; i++){
			int delta = Delta[i-1];
			for(int j=0; j<4; j++){
				table[X[j]][d] = F[i-1]*Y[(j+delta)%4];
				table[Y[(j+delta)%4]][d] = -F[i-1]*X[j];
			}
			d++;
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void fight(int a1, int b1, int a2, int b2, int a3, int b3, int f, int d, vector<vector<int> > &day)
{	
	day[a1][d]=f*b1, day[b1][d]=-f*a1;
	day[a2][d]=f*b2, day[b2][d]=-f*a2;
	day[a3][d]=f*b3, day[b3][d]=-f*a3;
}
void type_m(vector<int> &a, vector<int> &b, int s, int f, vector<vector<int> > &day, int d)
{
	if(s==1) fight(a[1],b[1], a[2],b[2], a[3],b[3], f, d, day);
	else if(s==2) fight(a[1],b[2], a[2],b[3], a[3],b[1], f, d, day);
	else if(s==3) fight(a[1],b[3], a[2],b[1], a[3],b[2], f, d, day);
	else if(s==4) fight(a[1],b[1], a[2],b[3], a[3],b[2], f, d, day);
	else if(s==5) fight(a[1],b[3], a[2],b[2], a[3],b[1], f, d, day);
	else if(s==6) fight(a[1],b[2], a[2],b[1], a[3],b[3], f, d, day);
}
void normalSupergameTriangle(int above, int below, vector<int> &x, vector<int> &y, vector<vector<double> > &newmatrix, vector<vector<int> > &table, int d, int f)
{
	vector<int> a(4), b(4);
	
	int x1=1,y1=2,z1=3;
	
	a[1] = x[3*above-2]; 
	a[2] = x[3*above-1]; 
	a[3] = x[3*above];
	
	b[1] = y[3*below-2]; 
	b[2] = y[3*below-1]; 
	b[3] = y[3*below];
	
	double 
	s1 = newmatrix[a[1]][b[1]]+newmatrix[a[2]][b[2]]+newmatrix[a[3]][b[3]], 
	s2 = newmatrix[a[1]][b[2]]+newmatrix[a[2]][b[3]]+newmatrix[a[3]][b[1]],
	s3 = newmatrix[a[1]][b[3]]+newmatrix[a[2]][b[1]]+newmatrix[a[3]][b[2]],
	s4 = newmatrix[a[1]][b[1]]+newmatrix[a[2]][b[3]]+newmatrix[a[3]][b[2]],
	s5 = newmatrix[a[1]][b[3]]+newmatrix[a[2]][b[2]]+newmatrix[a[3]][b[1]],
	s6 = newmatrix[a[1]][b[2]]+newmatrix[a[2]][b[1]]+newmatrix[a[3]][b[3]];
		
	double temp1=s1+s2+s3-max(max(s1,s2),s3), temp2=s4+s5+s6-max(max(s4,s5),s6);
	
	if(temp1<=temp2){
		if(s1>=s2 && s1>=s3) x1=2,y1=1,z1=3;
		else if(s2>=s1 && s2>=s3) x1=1,y1=2,z1=3;
		else x1=1,y1=3,z1=2;
	}
	else{
		if(s4>=s5 && s4>=s6) x1=5,y1=4,z1=6;
		else if(s5>=s4 && s5>=s6) x1=4,y1=5,z1=6;
		else x1=4,y1=6,z1=5;
	}
	
	type_m(a,b,x1,f,table,d+0); type_m(a,b,x1,-f,table,d+3); 
	type_m(a,b,y1,f,table,d+1); type_m(a,b,y1,-f,table,d+4); 
	type_m(a,b,z1,f,table,d+2); type_m(a,b,z1,-f,table,d+5);
	
	vector<int>().swap(a), vector<int>().swap(b);
}
void leftSupergameTriangle(int n, int above, int below, vector<int> &x, vector<int> &y, vector<vector<int> > &table, int d, int f)
{
	//x[3*above-2], x[3*above-2], x[3*above-1], x[n] ------- y[3*below-2], y[3*below-2], y[3*below-1], y[n]
	int X[4]={x[3*above-2], x[3*above-1], x[3*above], x[n]};
	int Y[4]={y[3*below-2], y[3*below-1], y[3*below], y[n]};
	
	for(int i=1; i<=3; i++){
		for(int j=0; j<4; j++){
			table[X[j]][d] = f*Y[(j+i)%4];
			table[Y[(j+i)%4]][d] = -f*X[j];
		}
		d++;
	}
	f = -f;
	for(int i=1; i<=3; i++){
		for(int j=0; j<4; j++){
			table[X[j]][d] = f*Y[(j+i)%4];
			table[Y[(j+i)%4]][d] = -f*X[j];
		}
		d++;
	}
}
void lastSupergame(int n, vector<int> &x, vector<int> &y, vector<vector<int> > &table, int f)
{
	int m=(n-1)/3;
	
	table[x[n]][1] = -f*y[n];
	table[y[n]][1] = f*x[n];
	
	for(int i=1; i<=m; i++){
		int above = i, below = 2*i-1;
		while(below>m) below -= m;
		
		int X[3]={x[3*above-2], x[3*above-1], x[3*above]};
		int Y[3]={y[3*below-2], y[3*below-1], y[3*below]};
		
		for(int i=0; i<3; i++){
			table[X[i]][1] = -f*Y[i];
			table[Y[i]][1] = f*X[i];
		}
	}
	
	table[x[n]][2*n] = f*y[n];
	table[y[n]][2*n] = -f*x[n];
	
	for(int i=1; i<=m; i++){
		int above = i, below = 2*i-1;
		while(below>m) below -= m;
		
		int X[3]={x[3*above-2], x[3*above-1], x[3*above]};
		int Y[3]={y[3*below-2], y[3*below-1], y[3*below]};
		
		for(int i=0; i<3; i++){
			table[X[i]][2*n] = f*Y[i];
			table[Y[i]][2*n] = -f*X[i];
		}
	}

}
double construction(int n, vector<int> &x, vector<int> &y, vector<vector<double> > &newmatrix, vector<vector<int> > &table, int direction)
{
	int m=(n-1)/3;

	for(int i=1; i<=m; i++)
	{
		//printf("Time slot %d: ", i);
		for(int j=1; j<=m; j++)
		{
			// above plays a super-team with below
			int above = i+j-1, below = 2*i+j-2;
			while(above>m) above -= m;
			while(below>m) below -= m;
			
			if(j==1)
			{
				leftSupergameTriangle(n, above, below, x, y, table, 6*i-4, direction);
			}
			else
			{
				normalSupergameTriangle(above, below, x, y, newmatrix, table, 6*i-4, direction);
			}
		}
	}
	lastSupergame(n, x, y, table, direction);
	
	//show_table(n, table);
	check(n, table);
	double ans = dist(n, table, newmatrix);
	
	//printf("%lf\n", ans);
	//scanf("%d", &m);
	return ans;	
}
void swap_team(int a, int b, int n, vector<int> &x, vector<int> &y)
{
	if(a<n){
		int tem = x[a];
		x[a] = x[b];
		x[b] = tem;
	}
	else{
		a -= n; b -= n;
		int tem = y[a];
		y[a] = y[b];
		y[b] = tem;
	}
	
}
int swap_normal_team(int n, vector<int> &x, vector<int> &y, vector<vector<double> > &newmatrix, vector<vector<int> > &currentTable, vector<vector<int> > &bestTable, int direction)
{
	int f_improvement=0;
	double minResult = construction(n, x, y, newmatrix, bestTable, direction); printf("%lf, ", minResult);
	vector<pair<int, int> > list_n;
	for(int i=1; i<=n; i++) {
		for(int j=i+1; j<=n; j++) {
			list_n.push_back(make_pair(i,j));
			list_n.push_back(make_pair(i+n,j+n));
		}
	}
	
	int flag=1;
	
	while(flag)
	{
		random_shuffle(list_n.begin(), list_n.end());
		
		for(int i=0; i<list_n.size(); i++) 
		{
			swap_team(list_n[i].first, list_n[i].second, n, x, y);
			
			double temResult = construction(n, x, y, newmatrix, currentTable, direction);
			
			if(temResult<minResult-0.001) 
			{
				minResult = temResult, flag = 2; 
				for(int i=1; i<=2*n; i++){
					for(int j=1; j<=2*n; j++){
						bestTable[i][j] = currentTable[i][j];
					}
				}
				printf("%lf, ", minResult);
				f_improvement=1;
			} 
			else swap_team(list_n[i].first, list_n[i].second, n, x, y);
		}
		
		flag--;
	}
	
	printf("\n");
	
	vector<pair<int, int> >().swap(list_n);
	
	return f_improvement;
}

int main() 
{
	//srand(time(0));
	clock_t s_clock, t_clock; 
	s_clock = clock();
    // calculate the distance of home venues
    int n = 16;
    vector<vector<double> > coordinates(2*n+1, vector<double>(2+1));
    FILE *fp = NULL;
    fp = fopen("locations.txt","r");
    
	for(int i=1; i<=2*n; i++){
		fscanf(fp,"%lf %lf", &coordinates[i][1], &coordinates[i][2]);
	}

	vector<vector<double> > newmatrix(2*n+1, vector<double>(2*n+1));
	for(int i=1; i<=2*n; i++){
		for(int j=1; j<=2*n; j++){
			newmatrix[i][j] = haversineDistance(coordinates[i][1], coordinates[i][2], coordinates[j][1], coordinates[j][2]);
		}
	}
	// find x_16 and y_16
	vector<double> costx(n+1), costy(n+1);
	for(int i=1; i<=n; i++){
		for(int j=1; j<=n; j++){
			costx[i] += newmatrix[i][j+n];
			costy[i] += newmatrix[n+i][j];
		}
	}
	
	double mincostx = 9999999.0, mincosty = 99999999.0;
	int x16 = 0, y16 = 0;	
	for(int i=1; i<=n; i++)
	{
		if(costx[i]<mincostx){
			mincostx = costx[i];
			x16 = i;
		}
		
		if(costy[i]<mincosty){
			mincosty = costy[i];
			y16 = i + n;
		}
	}
	
	vector<int> temx, temy, x(n+1), y(n+1);
	for(int i=1; i<=n; i++){
		if(i != x16) temx.push_back(i);
		if(i+n != y16) temy.push_back(i+n);
	}
	temx.push_back(x16);
	temy.push_back(y16);
	
	for(int i=1; i<=n; i++){
		x[i] = temx[i-1];
		y[i] = temy[i-1];
	}
	// find the best solution
	int direction = 1;
	vector<vector<int> > currentTable(2*n+1, vector<int>(2*n+1)), bestTable(2*n+1, vector<int>(2*n+1));
	swap_normal_team(n, x, y, newmatrix, currentTable, bestTable, direction);
	
    printf("After local search, the best table is:\n");
	show_table(n, bestTable);
 	printf("The total traveling distance is:%lf\n", dist(n, bestTable, newmatrix));
 	
	t_clock = clock();
 	printf("The running time is:%lf\n", double(t_clock-s_clock)/CLOCKS_PER_SEC);	
 	printf("Click to exist");
 	getchar(); 
    return 0;
}
