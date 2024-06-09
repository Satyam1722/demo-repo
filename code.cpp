#include<bits/stdc++.h>
// #include<ext/pb_ds/assoc_container.hpp>
// #include<ext/pb_ds/tree_policy.hpp>

using namespace std;
// using namespace __gnu_pbds;

// typedef tree<int, null_type, less<int>, rb_tree_tag, tree_order_statistics_node_update> oset; 
// find_by_order, order_of_key
typedef long long int ll;
typedef long double ld;
#define ff first
#define ss second
#define pb push_back

void _print(int a) {cerr<<a;}
void _print(ll a) {cerr<<a;}
void _print(char a) {cerr<<a;}
void _print(string a) {cerr<<a;}
void _print(double a) {cerr<<a;}

// Function  Declartion
void solve(void);

template<class T> void _print(vector<T> v){ cerr<<"[ "; for(T x : v){ _print(x); cerr<<" ";} cerr<<"]";}
template<class T , class R> void _print(pair<T , R> v){ cerr<<"[ "; cerr<<v.ff<<" "<<v.ss<<" "; cerr<<"]";}
template<class T> void _print(set<T> v){ cerr<<"[ "; for(T x : v) cerr<<x<<" "; cerr<<"]";}
template<class T> void _print(unordered_set<T> v){ cerr<<"[ "; for(T x : v) cerr<<x<<" "; cerr<<"]";}
template<class T> void _print(multiset<T> v){ cerr<<"[ "; for(T x : v) cerr<<x<<" "; cerr<<"]";}
template<class T> void _print(T a[] , int size){cerr<<"[ "; for(int i = 0 ; i < size ;++i) cerr<<a[i]<<" "; cerr<<"]"; }
template<class T> void _print(stack<T> v){ cerr<<"[ "; while(v.size() > 0) { cerr<<v.top()<<" "; v.pop() ; } cerr<<"]";}
template<class T , class R> void _print(map<T,R> m){cerr<<"[ "; for(auto x : m) cerr<<x.first<<" "<<x.second<<endl; cerr<<"]";}
template<class T , class R> void _print(vector<pair<T,R>> m){cerr<<"[ "; for(auto x : m) cerr<<x.first<<" "<<x.second<<endl; cerr<<"]";}


#define readvll(x , a) vector<ll> x(a); for(int i = 0 ; i < a ; i++) cin >> x[i];
#define vll vector<ll>
#define vi vector<int>
#define rpq priority_queue<ll,vector<ll>,greater<ll>> 
#define pq priority_queue<ll>
#define vllp vector<pair<ll,ll>>
#define full(a) a.begin(),a.end()
#define sor(a) sort(full(a))
#define rev(a) a.rbegin(),a.rend()
#define pb push_back
#define dg cout<<"hello world"<<endl
#define dg1 cout<<"hello world1"<<endl
#define dg2 cout<<"hello world2"<<endl
#define endl "\n"

#ifndef ONLINE_JUDGE
#define debug(x) cerr<<#x<<" "; _print(x); cerr<<endl;
#else
#define debug(x)
#endif

ll min(vll v){	ll mn = 1e18 ; for( ll i : v ) {  if( mn > i) mn=i; } return mn;}
ll max(vll v){	ll mn = -1e18 ; for( ll i : v ) {  if( mn < i) mn=i; } return mn;}
int dx[]={-1,0,1,0};
int dy[]={0,1,0,-1};
int knightx[]={-1,-2,-2,-1,1,2,2,1};
int knighty[]={-2,-1,1,2,2,1,-1,-2};

vector<int> isPrime;
void  sieve(int n){
	n++;
	vector<int> v(n , 1);
	v[0] = v[1] = 0;
	for(int i = 2 ; i <= sqrt(n) ; i++){
		if(v[i]){
			for(int j = 2 * i ; j < n ; j += i){
				v[j] = 0;
			}
		}
	}
	isPrime = v;
}

// fact , inverese factorial and inverse modulo
vector<ll> invMod, invFact, fact;
const ll mod = 1e9+7;
const ll mod1 = 998244353;
ll lcm(ll a , ll b){	return a*b/__gcd(a,b);		}
ll mul(ll a, ll b){    return ((a % mod)*(b % mod)) % mod;}
ll add(ll a, ll b){    return ((a % mod)+(b % mod)) % mod;}
ll power(ll x, ll y)
{
	long long int temp;
	if (y == 0)
		return 1;
	temp = power(x, y / 2);
	if (y % 2 == 0)
		return mul(temp, temp);
	else
		return mul(x, mul(temp, temp));
}



void invModf(int n){
	invMod = vector<ll>(n + 1);
	invMod[1]  = 1;
	for(int i = 2; i <=n; i++){
		invMod[i] = (mod - mul((mod / i), invMod[mod % i]))%mod;
	}
}
void invFactf(int n){
	invFact = vector<ll>(n + 1);
	invFact[1] = invFact[0] = 1;
	for(int i = 2; i <= n; i++){
		invFact[i] = mul(invFact[i - 1], invMod[i]);
	}
}
void factf(int n){
	fact = vector<ll>(n + 1);
	fact[1] = fact[0] = 1;
	for(int i = 2; i <= n; i++){
		fact[i] = mul(i , fact[i - 1]);
	}
}
int nCr(int n, int r){
	return mul(fact[n], mul(invFact[r], invFact[n- r]));
}
void initNCR(int n){
	invModf(n);
	invFactf(n);
	factf(n);
}
int nDigit(ll num , int base){
	int dig = 0;
	while(num > 0){
		num /= base;
		dig++;
	}
	return dig;
}
ll ETF(ll n) {  // gives the count of positive integers less than or equal to n that are coprime
	ll result = n;  // Initialize result with n

	// Check for each prime factor
	for (int i = 2; i * i <= n; ++i) {
		if (n % i == 0) {
			while (n % i == 0) {
				n /= i;
			}
			result -= result / i;
		}
	}

	// If n is a prime number greater than 1
	if (n > 1) {
		result -= result / n;
	}

	return result;
}
bool getbit(long long int n, int i) {
	return n & (1LL << (i));
}
long long int setbit(long long int n, int i) {
	return (n | (1LL << (i)));
}
long long int clearbit(long long int n, int i) {
	return n & (~(1LL << (i)));
}
long long int togglebit(long long int n, int i) {
	return n ^ (1LL << (i));
}
int countbit(long long int n){
	return __builtin_popcount(n);
}

map<ll,ll> primeFac(ll n)
{
	map<ll ,ll> v;
	while (n % 2 == 0)
	{
		v[2]++;
		n = n/2;
	}
	ll sqrtNum = static_cast<ll>(sqrt(n));
	for (ll i = 3; i <= sqrtNum; i = i + 2)
	{
		while (n % i == 0)
		{
			v[i]++;
			n = n/i;
		}
	}
	if (n > 2)
	v[n]++;

	return v;
}

vector<ll>divisor(ll number) {
	vector<ll> divisors;
	ll sqrtNumber = static_cast<ll>(sqrt(number));
	
	for (int i = 1; i <= sqrtNumber; ++i) {
		if (number % i == 0) {
			divisors.push_back(i);
			if (i != number / i) {
				divisors.push_back(number / i);
			}
		}
	}
	
	return divisors;
}

bool isprime(ll n)
{

	ll sqrtNum = static_cast<ll>(sqrt(n));
	for (ll i = 2; i <= sqrtNum; i++)
	{
		if(n % i == 0)
		{
			return false;
		}
	}

	return true;
}

vector<ll> _div;
void nDivisor(int n){
	n++;
	vector<ll> v(n , 0);
	v[0] = v[1] = 0;
	for(int i = 2 ; i < n; i++){
		if((i & (i - 1)) == 0){
			for(int j = i ; j < n ; j += i){
				v[j]++;
			}
		}
	}
	_div = v;
}
// Function to calculate the modular inverse of 'a' modulo 'm'
ll modInverse(ll a) {
	return power(a,mod-2);
}
ll revn(ll num) {
	ll reversedNum = 0;

	while (num != 0) {
		ll digit = num % 10;
		reversedNum = reversedNum * 10 + digit;
		num /= 10;
	}

	return reversedNum;
}

vector<ll> bina(ll n , int count){
	vector<ll> v(count,0);
	for(int i=0; n>0; i++)    
	{    
		v[i] = (n%2);    
		n= n/2;  
	}
	return v;
}
ll rbin(vector<ll> &v){
	ll ans = 0;
	int n = v.size();
	for(int i = 0 ; i < n; i++){
		ans += power(2,i)*v[i];
	}
	return ans;
}

vector<int> computeLPS(const string& pattern) {
	int m = pattern.length();
	vector<int> lps(m, 0);
	int j = 0;

	for (int i = 1; i < m; ++i) {
		while (j > 0 && pattern[i] != pattern[j]) {
			j = lps[j - 1];
		}
		if (pattern[i] == pattern[j]) {
			lps[i] = ++j;
		}
	}
	return lps;
}

int KMP(const string& text, const string& pattern) {
	int n = text.length();
	int m = pattern.length();
	vector<int> lps = computeLPS(pattern);
	int j = 0;

	for (int i = 0; i < n; ++i) {
		while (j > 0 && text[i] != pattern[j]) {
			j = lps[j - 1];
		}
		if (text[i] == pattern[j]) {
			if (j == m - 1) {
				// cout << "Pattern found at index " << i - m + 1 << endl;
				return i-m+1;
				j = lps[j];
			} else {
				++j;
			}
		}
	}

	return -1;
}

vector<int> computeZ(const string& str) {
	int n = str.length();
	vector<int> Z(n, 0);
	int left = 0, right = 0;

	for (int i = 1; i < n; ++i) {
		if (i <= right) {
			Z[i] = min(Z[i - left], right - i + 1);
		}
		while (i + Z[i] < n && str[Z[i]] == str[i + Z[i]]) {
			++Z[i];
		}
		if (i + Z[i] - 1 > right) {
			left = i;
			right = i + Z[i] - 1;
		}
	}
	return Z;
}

vector<int> ZAlgorithm(const string& text, const string& pattern) {
	string concat = pattern + "$" + text;
	vector<int> Z = computeZ(concat);
	int patternLength = pattern.length();
	vector<int> indices;

	for (int i = patternLength + 1; i < concat.length(); ++i) {
		if (Z[i] == patternLength) {
			indices.push_back(i - patternLength - 1);
		}
	}
	return indices;
}


bool cmp(pair<ll ,ll>& a,pair<ll ,ll>& b)
{
	if(a.ss<b.ss)
	return true;	
	else if(a.ss==b.ss && a.ff>b.ff)
	return true;
	else
	return false;
}


vector<string> tokens;	

void tokenize(string s,char delimiter)
{

	istringstream iss(s);
	string token;

	while (getline(iss, token, delimiter)) {
		tokens.push_back(token);
	}
}

std::vector<ll> split(const std::string& input, char delimiter) {
	std::vector<ll> tokens;
	std::istringstream stream(input);
	std::string token;

	while (std::getline(stream, token, delimiter)) {
		tokens.push_back(stoll(token));
	}

	return tokens;
}	

struct Node{
	Node* links[2];
	int count=0;

	bool isContain(int bit)
	{
		return links[bit]!=NULL;
	}

	void put(int bit,Node *node)
	{
		links[bit]=node;
	}

	Node* get(int bit)
	{
		return links[bit];
	}

void setEnd()
{
		count++;
}

void unsetEnd()
{
		count--;
}

int get()
{
		return count;
}
};

class Trie{

	Node *root;
	public:

	Trie(){
		// Write your code here.
		root=new Node();
	}

	void insert(int num){
		// Write your code here.
		Node *node=root;

		for(int i=31;i>=0;i--)
		{
			int bit=(num>>i)&1;

			if(!node->isContain(bit))
			{
				node->put(bit,new Node());
			}
			node=node->get(bit);
		}
		node->setEnd();
	}

	int maxXor(int num){
		// Write your code here.
		Node* node=root;
		int ans=0;

		for(int i=31;i>=0;i--)
		{
			int bit=(num>>i)&1;

			if(node->isContain(1-bit))
			{
				ans=ans|(1<<i);
				node=node->get(1-bit);
			}
			else
			{
				node=node->get(bit);
			}
		}

		return ans;
	}

	int minXor(int num)
	{
		Node *node=root;
		int ans=0;

		for(int i=31;i>=0;i--)
		{
			int bit=(num>>i)&1;

			if(!node->isContain(bit))
			{
				ans=ans|(1<<i);
				node=node->get(1-bit);
			}
			else
			{
				node=node->get(bit);
			}
		}

		return ans;
	}


	void erase(int num){
		// Write your code here.
		Node* node=root;
		stack<pair<Node*,int>> st;

		for(int i=31;i>=0;i--)
		{
			int bit=(num>>i)&1;
			st.push({node,bit});
			node=node->get(bit);
		}
		node->unsetEnd();

		if(node->get()!=0)
		return;

		while(!st.empty())
		{
			auto it=st.top();
			st.pop();

			Node *parent=it.first;
			int bit=it.second;
			Node *child=it.first->get(bit);

			if(child->isContain(1) || child->isContain(0) || child->get()!=0)
			break;

			it.first->put(bit,NULL);
			delete child;
		}
	}
	};

int main()
{
	#ifndef ONLINE_JUDGE
		freopen("input.txt","r",stdin);
		freopen("output.txt","w",stdout);
		freopen("error.txt","w",stderr);
	#endif

	ios_base::sync_with_stdio(false);cin.tie(NULL);cout.tie(NULL);

	// sieve(1000001);
	// ndivisor(200010);
	// initNCR(200001);
	// calc();

	// ll t;
	// cin>>t;
	// while(t--)
	{
		solve();
	}

	return 0;
}

void dfs(int u,int p,vector<vi> &adj,vector<vi> &parent,vi &level)
{
	parent[0][u]=p;

	if(p!=-1)
	{
		level[u]=level[p]+1;
	}

	for(auto v:adj[u])
	{
		if(v!=p)
		{
			dfs(v,u,adj,parent,level);
		}
	}
}

void dfs(int u,int p,vector<vi> &adj,vll &sum,vll &dist,vll &a)
{
	for(auto v:adj[u])
	{
		if(v!=p)
		{
			dfs(v,u,adj,sum,dist,a);
			sum[u]+=sum[v]+a[v];
			dist[u]+=dist[v]+sum[v]+a[v];
		}
	}
}

void reroot(int u,int v,vector<vi> &adj,vi &score)
{
	score[u]-=max(0,score[v]);
	score[v]+=max(0,score[u]);
}

void dfs(int u,int p,vector<vi> &adj,vi &score,vi &ans)
{
	ans[u]=score[u];

	for(auto v:adj[u])
	{
		if(v!=p)
		{
			reroot(u,v,adj,score);
			dfs(v,u,adj,score,ans);
			reroot(v,u,adj,score);
		}
	}
}

int findParent(int u,int k,vector<vi> &parent,int maxPow)
{
	int ans=u;
	// cout<<ans<<endl;

	for(int i=0;i<maxPow;i++)
	{
		if(k&(1<<i))
		{
			ans=parent[i][ans];
			// cout<<ans<<" "<<i<<endl;
			if(ans==-1)
			return ans;
		}
	}

	return ans;
}



int lca(int a,int b,vector<vi> &parent,vi &level,int maxPow)
{
	if(level[a]>level[b])
	{
		a=findParent(a,level[a]-level[b],parent,maxPow);
	}
	else
	{
		b=findParent(b,level[b]-level[a],parent,maxPow);
	}

	if(a==b)
	return a;

	int low=0;
	int high=level.size()-1;

	int ans=-1;

	while(low<=high)
	{
		int mid=low+(high-low)/2;

		int u=findParent(a,mid,parent,maxPow);
		int v=findParent(b,mid,parent,maxPow);

		if(u==v)
		{
			ans=u;
			high=mid-1;
		}
		else
		{
			low=mid+1;
		}

	}

	return ans;
}

int lca2(int a,int b,vector<vi> &parent,vi &level,int maxPow)
{
	if(level[a]>level[b])
	{
		a=findParent(a,level[a]-level[b],parent,maxPow);
	}
	else
	{
		b=findParent(b,level[b]-level[a],parent,maxPow);
	}

	if(a==b)
	return a;

	for(int i=maxPow-1;i>=0;i--)
	{
		int u=parent[i][a];
		int v=parent[i][b];

		if(u==v) continue;

		a=u;
		b=v;
	}

	return parent[0][a];
}


int getdist(int u,int k,int maxPow,vector<vi> &parent,vector<vi> &dist)
{
	int ans=0; 

	for(int i=maxPow-1;i>=0;i--)
	{
		if(k&(1<<i))
		{
			ans+=dist[i][u];
			u=parent[i][u];
		}
	}

	return ans;
}


void dfs(int u,int p,vector<vi> &adj,vi &arr,vi &intime,vi &outime,int &time)
{
	intime[u]=time++;
	arr.pb(0);

	for(auto v:adj[u])
	{
		if(v!=p)
		{
			dfs(v,u,adj,arr,intime,outime,time);
		}
	}
	outime[u]=time++;
	arr.pb(0);
}


void helper(ll sr,vector<vector<pair<ll,ll>>> &adj,vll &dist)
{
	dist[sr]=0;
	set<pair<ll,ll>> st;
	st.insert({dist[sr],sr});

	while(!st.empty())
	{
		pair<ll,ll> node=*st.begin();
		st.erase(st.begin());

		ll d=node.ff;
		ll u=node.ss;

		for(auto it:adj[u])
		{
			ll v=it.ff;
			ll w=it.ss;

			if(d+w<dist[v])
			{
				st.erase({dist[v],v});
				dist[v]=d+w;
				st.insert({dist[v],v});
			}
		}
	}

}

void helper(ll u,vector<vll> &adj,vll &dist)
{
	int n=adj.size();
	dist[u]=0;
	set<pair<ll,ll>> st;
	st.insert({dist[u],u});

	while(!st.empty())
	{
		pair<ll,ll> node=*st.begin();
		st.erase(st.begin());

		ll d=node.ff;
		ll u=node.ss;

		for(int v=0;v<n;v++)
		{
			ll w=adj[u][v];

			if(d+w<dist[v])
			{
				st.erase({dist[v],v});
				dist[v]=d+w;
				st.insert({dist[v],v});
			}
		}
	}

}

void solve()
{
	ll n;
	cin>>n;
	
	vector<vll> adj;

	for(int i=0;i<n;i++)
	{
		readvll(v,n);
		adj.pb(v);
	}

	vector<vll> dist(n,vll(n,INT_MAX));

	readvll(v,n);
	reverse(full(v));

	vll ans;

	vll removed(n,true);

	
	for(auto it:v)
	{ 
		int k=--it;
		removed[k]=false;

		for(int i=0;i<n;i++)
		{
			dist[i][k]=min(dist[i][k],adj[i][k]);
			dist[k][i]=min(dist[k][i],adj[k][i]);
		}

		ll sum=0;
		for(int i=0;i<n;i++)
		{
			for(int j=0;j<n;j++)
			{
				dist[i][j]=min(dist[i][j],dist[i][k]+dist[k][j]);
				// cout<<i<<" "<<j<<endl;
				// cout<<dist[i][j]<<endl;

				if(removed[i] || removed[j] || i==j)
				continue;
				
				// cout<<i<<" "<<j<<endl;
				// cout<<removed[i]<<" "<<removed[j]<<endl;
				// cout<<dist[i][j]<<endl;
				sum+=dist[i][j];
			}
		}

		// cout<<sum<<endl;

		ans.pb(sum);
	}

	reverse(full(ans));

	for(auto it:ans)
	cout<<it<<" ";
	cout<<endl;

}


