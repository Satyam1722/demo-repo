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


void dfs(int u,int p,vector<vi> &adj,vi &score)
{
	for(auto v:adj[u])
	{
		if(v!=p)
		{
			dfs(v,u,adj,score);
			score[u]+=max(0,score[v]);
		}
	}
}

void solve()
{
	ll a,b;
	cin>>a>>b;

	a=a-b;
	b=a+2*b;
