/*
  Experiment for graph visalization in 2D
*/
int N;
int[][] edge;
PVector[] p;
PVector[] dp;
float[][] norm;

float[][] ze;
float[][] zd;
float[][] dEdd;

void settings()
{
  size(640, 640);
}

void setup()
{
  colorMode(HSB, 1);
  N = 16;
  create_random_graph();
}

void draw()
{
  background(0);
  draw_graph();
  update();
}

void create_random_graph()
{
  edge = new int[N][N];
  p = new PVector[N];
  dp = new PVector[N];
  norm = new float[N][N];
  ze = new float[N][N];
  zd = new float[N][N];
  dEdd = new float[N][N];
  for(int n = 0; n < N; ++n)
  {
    p[n] = new PVector(random(width/4, width*3/4), random(height/4, height*3/4));
    dp[n] = new PVector();
  }
  for(int n = 0; n < N; ++n)
    for(int m = n + 1; m < N; ++m)
    {
      if(random(1) > 0.5)
        edge[n][m] = (int)constrain(randomGaussian() * 50 + 50, 1, 100);
      norm[n][m] = dist(p[n].x, p[n].y, p[m].x, p[m].y);
    }
}

// this assume that distance will be synchronized to
// edge weight. also this assumes edge weight >= 0.
void update()
{
  int cnt = 0;
  
  // refresh distance info
  for(int n = 0; n < N; ++n)
    for(int m = n + 1; m < N; ++m)
      norm[n][m] = dist(p[n].x, p[n].y, p[m].x, p[m].y);
  
  // compute the summation of weights, distance
  float edge_sum = 0;
  float norm_sum = 0;
  for(int n = 0; n < N; ++n)
    for(int m = n + 1; m < N; ++m)
      if(edge[n][m] > 0)
      {
        edge_sum += edge[n][m];
        norm_sum += norm[n][m];
      }
    
  // normalize
  // trace error
  float err = 0;
  for(int n = 0; n < N; ++n)
    for(int m = n + 1; m < N; ++m)
    {
      if(edge[n][m] > 0)
      {
        ++cnt;
        ze[n][m] = edge[n][m] / edge_sum;
        zd[n][m] = norm[n][m] / norm_sum;
        err += sq(zd[n][m] - ze[n][m]);
      }
    }
  println(err);
  
  // compute gradient of norm_sum = dEdS
  float dEdS = 1e-10 * (norm_sum - cnt * 200);
  for(int n = 0; n < N; ++n)
    for(int m = n + 1; m < N; ++m)
      if(edge[n][m] > 0)
        dEdS -= 2 * (zd[n][m] - ze[n][m]) * zd[n][m] / norm_sum;
      
  // compute dE/dd (each edge's distance's gradient)
  for(int n = 0; n < N; ++n)
    for(int m = n + 1; m < N; ++m)
      if(edge[n][m] > 0)
        dEdd[n][m] = 2 * (zd[n][m] - ze[n][m]) / norm_sum + dEdS;
  
  // compute dE/dx, dE/dy of each vertices
  for(int n = 0; n < N; ++n)
  {
    dp[n].x = 0;
    dp[n].y = 0;
    for(int m = 0; m < N; ++m)
    {
      if(m > n)
      {
        if(edge[n][m] == 0)
          continue;
        dp[n].x += dEdd[n][m] * (p[n].x - p[m].x) / (norm[n][m] + 1e-3);
        dp[n].y += dEdd[n][m] * (p[n].y - p[m].y) / (norm[n][m] + 1e-3);
      }
      else if(m < n)
      {
        if(edge[m][n] == 0)
          continue;
        dp[n].x += dEdd[m][n] * (p[n].x - p[m].x) / (norm[m][n] + 1e-3);
        dp[n].y += dEdd[m][n] * (p[n].y - p[m].y) / (norm[m][n] + 1e-3);
      }
    }
  }
  
  // update gd
  // you can change this to another adaptive algorithm
  for(int n = 0; n < N; ++n)
  {
    p[n].x -= dp[n].x * cnt * 20000;
    p[n].y -= dp[n].y * cnt * 20000;
  }
}

void draw_graph()
{
  for(int n = 0; n < N; ++n)
    for(int m = n + 1; m < N; ++m)
      if(edge[n][m] > 0)
      {
        stroke(0.1, edge[n][m] / 100.0, 1);
        line(p[n].x, p[n].y, p[m].x, p[m].y);
      }
  fill(0, 1, 0.95);
  noStroke();
  for(int n = 0; n < N; ++n)
    ellipse(p[n].x, p[n].y, 10, 10);
}

void keyReleased()
{
  if(key == ' ')
    create_random_graph();
}
