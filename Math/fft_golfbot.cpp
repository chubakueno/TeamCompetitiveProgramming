#include <bits/stdc++.h>
using namespace std;
typedef complex<double> cld;
cld wlen_pw[262144];
vector<cld> v1;
vector<int> vec;
int rev[2*262144];
//invertir el orden de los bits
void calc_rev(int n) {
	int pw=1<<(31-__builtin_clz(n-1));
	rev[0]=0;
	rev[1]=pw;
	for (int i=2; i<n; i+=2)
		rev[i+1]=(rev[i] = rev[i>>1]>>1)|pw;
}
void fft (cld a[], int n, bool invert){
	for (int i=0; i<n; ++i)
		if (i < rev[i])
			swap (a[i], a[rev[i]]);
	for (int len=2; len<=n; len<<=1) {
		double ang = 2*M_PI/len * (invert?-1:+1);
		int len2 = len>>1;
		cld wlen (cos(ang), sin(ang));
		wlen_pw[0] = 1;
		for (int i=1; i<len2; ++i)
			wlen_pw[i] = wlen_pw[i-1] * wlen;
		for (int i=0; i<n; i+=len) {
			cld t,
				*pu = a+i,
				*pv = a+i+len2, 
				*pu_end = a+i+len2,
				*pw = wlen_pw;
			for (; pu!=pu_end; ++pu, ++pv, ++pw) {
				t = *pv * *pw;
				*pv = *pu - t;
				*pu += t;
			}
		}
	}
	if (invert)
		for (int i=0; i<n; ++i)
			a[i] /= n;
}
int main(){
	int n;
	while(scanf("%d",&n)==1){
		int mv=0;
		vec.clear();
		for(int i=0;i<n;++i){
			int x;
			scanf("%d",&x);
			vec.push_back(x);
			mv=max(mv,x);
		}
		int sz = 1;
		while (sz < mv+1)  sz <<= 1;
		sz<<=1;
		calc_rev(sz);
		v1.assign(sz,0);
		v1[0]=1;
		for(int i=0;i<vec.size();++i)
			v1[vec[i]]=1;
		fft (&v1[0],sz, false);
		for (int i=0; i<sz; ++i)
			v1[i] *= v1[i];
		fft (&v1[0],sz, true);
		int m;
		scanf("%d",&m);
		int ans=0;
		for(int i=0;i<m;++i){
			int x;
			scanf("%d",&x);
			if(x<sz&&(int)(v1[x].real() + 0.5)) ++ans;
		}
		printf("%d\n",ans);
	}
}