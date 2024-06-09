#include<bits/stdc++.h>
using namespace std;
#define int long long

int n,k;
int arr[100100]; 

bool check(int x)
{
    int no_of_points = 0;
    for(int i=1; i<n ; i++)
    {
        no_of_points += (((arr[i]-arr[i-1])+x-1)/x) - 1;
	}
    return no_of_points <= k;
}
    
void solve()
{
    cin>>n>>k;
    for(int i=0; i<n; i++)
    {
        cin>>arr[i];
    }
    sort(arr, arr+n);
    int lo = 0;
    int hi = arr[n-1]-arr[0];
    int ans = -1;
    
    while(lo<=hi)
    {
        int mid = lo+(hi-lo)/2;
        if(check(mid))
        {
            ans = mid;
            hi= mid-1;
        }
        else
        {
            lo = mid+1;
        }
    }
    cout << ans << endl;
}

signed main()
{
    #ifndef ONLINE_JUDGE
		freopen("input.txt","r",stdin);
		freopen("output.txt","w",stdout);
		freopen("error.txt","w",stderr);
	#endif

    int t;
    cin>>t;
    
    while(t--)
        solve();

}