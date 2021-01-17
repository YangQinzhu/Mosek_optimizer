////
//  Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
//
//  File:      facility_location.cc
//
//  Purpose: Demonstrates a small one-facility location problem.
//
//  Given 10 customers placed in a grid we wish to place a facility
//  somewhere so that the total sum of distances to customers is
//  minimized.
//
//  The problem is formulated as a conic optimization problem as follows.
//  Let f=(fx,fy) be the (unknown) location of the facility, and let
//  c_i=(cx_i,cy_i) be the (known) customer locations; then we wish to
//  minimize
//      sum_i || f - c_i ||
//  where
//      ||.||
//  denotes the euclidian norm.
//  This is formulated as
//
//  minimize   sum(d_i)
//  such that  d_i ^ 2 > tx_i ^ 2 + ty_i ^ 2, for all i
//             tx_i = cx_i - fx, for all i
//             ty_i = cy_i - fy, for all i
//             d_i > 0, for all i
////

#include <memory>
#include <iostream>
#include <iomanip>
#include "fusion.h"
#include <cmath>

using namespace mosek::fusion;
using namespace monty;

int main(int argc, char ** argv)
{
  auto facilityloc = new_array_ptr<double, 2>(
  { 
    { 5., 11.},
    {12.,  8.},
    { 9., 11.},
    { 5., 15.},
    { 0., 16.},
    { 1., 12.},
    { 7., 13.},
    { 6., 18.},
    { 5., 18.},
    {11., 10.}
  } );

  auto customerloc = new_array_ptr<double, 2>(
  { 
    {14., 18.},
    { 4.,  9.},
    {17.,  0.},
    {13.,  9.},
    { 9.,  7.},
    { 1.,  0.},
    {17.,  8.},
    {13., 19.},
    {15., 10.},
    { 8.,  7.},
    { 3.,  6.}

  } );

  auto facilityVolume = new_array_ptr<double, 2>(
    { {5},{1},{9}, {3}, {4}, {8}, {1}, {4}, {0}, {28} } 
    // { {5},{1},{9}, {3}, {4}, {8}, {1}, {4}, {0}, {27} } 
    // { {9}, {3}, {4}, {8}, {1}, {4}, {0}, {22} } //51
    // { {8}, {1}, {4}, {0}, {25} } //成功 38
    // { {8}, {1}, {4}, {3}, {8} }  //24 成功
    // {{10}, {5}, {8}, {1} } 
  );

  auto customerVolume = new_array_ptr<double, 2>(
    { {9},{2},{0}, {4}, {9}, {2}, {7}, {7}, {9}, {8}, {6} } //62
    // { {0}, {4}, {9}, {2}, {7}, {7}, {8}, {8}, {6} } // 51
    // {  {2}, {7}, {7}, {8}, {8}, {6} }  //成功38
    // { {3}, {1}, {7}, {3}, {4}, {6} }  //成功24
    // {{10}, {5}, {4}, {5}, } 
  );

  int Ncustomer = customerloc->size(0); //取出个数
  int Nfacility = facilityloc->size(0);


  auto distance = std::make_shared<ndarray<double,2>>(shape(Nfacility, Ncustomer));
  for(int s1=0; s1<Nfacility; s1++){
    for(int s2=0; s2<Ncustomer; s2++){
      //norm2
      double diffxpow = pow((*facilityloc)(s1,0)-(*customerloc)(s2,0),2);
      double diffypow = pow((*facilityloc)(s1,1)-(*customerloc)(s2,1),2);

      // std::cout << diffxpow << " " << diffypow << " " << sqrt(diffxpow+diffypow) << std::endl;
      (*distance)(s1,s2) = sqrt(diffxpow+diffypow);
    }
  }
        
  // auto distance = new_array_ptr<double, 2>(

  // {{8.0, 6.2},
  //   {2.3, 4.5} } 
  // );

  // std::cout << (*distance) << std::endl; //和python计算的结果先做对照
  std::cout<< "result of distance calculation:" <<std::endl;
  for(int ij=0; ij<Nfacility*Ncustomer; ij++){
    for(int j=1; j<Ncustomer; j++){
      if(ij == j*Ncustomer){ 
        std::cout << " " <<std::endl;
      }
    }
    std::cout << (*distance)[ij]<< " ";
  }
  std::cout<<std::endl; //加上空格美观

  // facilityVolume  直接输出是地址
  std::cout << "num of warehouse:" << Ncustomer << "\n" 
            << "num of factories:" << Nfacility << std::endl; //输出确认结果

  Model::t M = new Model("BetterTransition"); auto _M = finally([&]() { M->dispose(); });
  Variable::t tran = M->variable("transition", Set::make(Nfacility, Ncustomer), Domain::greaterThan(0.0));
  // Variable::t tran = M->variable("transition", Set::make(Nfacility, Ncustomer), Domain::inRotatedQCone());

  ///-----------------
  //i 个工厂 j个仓库

  //tran 的j列相加，得到的向量和 facilityVolume相等
  /*
  把限制名称去掉，不然会出现名字重复的报错，
  这里的slice +1 不会产生越界
  */

  for(int i=0; i<Nfacility; i++){
    M->constraint(Expr::sum(tran->slice(new_array_ptr<int,1>({i,0}), new_array_ptr<int,1>({i+1,Ncustomer}))), Domain::equalsTo((*facilityVolume)[i]));
  }

  for(int j=0; j<Ncustomer; j++){
    M->constraint(Expr::sum(tran->slice(new_array_ptr<int,1>({0,j}), new_array_ptr<int,1>({Nfacility,j+1}))), Domain::equalsTo((*customerVolume)[j]));
  }

  std::cout << "num of warehouse:" << Ncustomer << "\n";
  M->objective("total_cost", ObjectiveSense::Minimize, Expr::sum(Expr::dot(tran, distance)));
  // M->objective("total_cost", ObjectiveSense::Minimize, Expr::dot(tran, distance)); //点乘之后即求和
  std::cout << (*distance) << std::endl;
  // M->selectedSolution(SolutionType::Integer);
  M->solve();
  
  std::cout << "num of warehouse:" << Ncustomer << "\n";

  auto res = tran->level(); //C++11auto可以在声明变量的时候根据变量初始值的类型自动为此变量选择匹配的类型

  //输出最终结果
  std::cout<< "result of output:" <<std::endl;
  for(int ij=0; ij<Nfacility*Ncustomer; ij++){
    for(int j=1; j<Ncustomer; j++){
      if(ij == j*Ncustomer){ 
        std::cout << " " <<std::endl;
      }
    }
    std::cout << (*res)[ij]<< " ";
  }
  std::cout<<std::endl;
}

// std::cout << std::setprecision(2)
//           << "shape:" << (*res) << " " << std::endl;

// for(int ij=0; ij<Nfacility*Ncustomer; ij++){
//   for(int j=0; j<Ncustomer; j++){
//     if(ij == j*Ncustomer){ 
//       std::cout << "\n" << std::endl;
//     }
//     std::cout << res->slice(new_array_ptr<int,1>({0,j}), new_array_ptr<int,1>({Ncustomer,j+1}));
//   }
// }

// for(int j=0; j<Ncustomer; j++){
//   if(ij == j*Ncustomer){ 
//     std::cout << "\n" << std::endl;
//   }
//   std::cout << res->slice(new_array_ptr<int,1>({0,j}), new_array_ptr<int,1>({Ncustomer,j+1}));
// }

// int main(int argc, char ** argv)
// {
//   auto customerloc = new_array_ptr<double, 2>(
//   { {12.,  2. },
//     {15., 13. },
//     {10.,  8. },
//     { 0., 10. },
//     { 6., 13. },
//     { 5.,  8. },
//     {10., 12. },
//     { 4.,  6. },
//     { 5.,  2. },
//     { 1., 10. }
//   } );


//   int N = customerloc->size(0);  //客人数目
//   std::cout << N << *customerloc << std::endl; //输出确认结果
//   Model::t M = new Model("FacilityLocation"); auto _M = finally([&]() { M->dispose(); });
//   // Variable holding the facility location
//   Variable::t f = M->variable("facility", Set::make(1, 2), Domain::unbounded()); //工厂位置
//   // Variable defining the euclidian distances to each customer
//   Variable::t d = M->variable("dist", Set::make(N, 1), Domain::greaterThan(0.0)); //距离， 大于0
//   // Variable defining the x and y differences to each customer;
//   Variable::t t = M->variable("t", Set::make(N, 2), Domain::unbounded()); 
//   M->constraint("dist measure",
//                 Var::hstack(d, t), //维度拼接
//                 Domain::inQCone(N, 3));

//   Variable::t fxy = Var::repeat(f, N); //工厂位置的显示
//   M->constraint("xy diff", Expr::add(t, fxy), Domain::equalsTo(customerloc)); //t应该是 客人位置-工厂位置=t，所以求得此等式

//   M->objective("total_dist", ObjectiveSense::Minimize, Expr::sum(d)); //总距离求和最小

//   M->solve();

//   auto res = f->level(); //获取最终的结果
//   std::cout << std::setprecision(2)  //输出精度
//             << "Facility location = " << (*res)[0] << "," << (*res)[1] << std::endl;
// }

  //----------失败
  // M->constraint("facility volume", Expr::sum(tran, 1), Domain::equalsTo(1.0));//facilityVolume
  // M->constraint("customer volume", Expr::sum(tran, 0), Domain::equalsTo({{1.}})); //customerVolume
  // ---------------

  //--------尝试成功
  // M->constraint("customer volume", tran->index(new_array_ptr<int,1>({0,0})), Domain::equalsTo(1.)); //customerVolume
  // M->constraint("customer volume", Expr::sum(tran->slice(new_array_ptr<int,1>({0,0}), new_array_ptr<int,1>({1,2}))), Domain::equalsTo((*facilityVolume)[0])); //customerVolume
  //-----------