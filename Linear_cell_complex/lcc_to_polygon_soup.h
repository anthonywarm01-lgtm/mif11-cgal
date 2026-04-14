// Copyright (c) 2022 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of LCC-Lab.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
////////////////////////////////////////////////////////////////////////////////
#ifndef LCC_TO_POLYGON_SOUP_H
#define LCC_TO_POLYGON_SOUP_H

#include<array>
#include<queue>
#include<map>
#include<unordered_map>
#include<vector>
#include<CGAL/Kernel_traits.h>
#include<CGAL/Polygon_mesh_processing/autorefinement.h>
#include<CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include<CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include<CGAL/Polygon_mesh_processing/internal/Corefinement/predicates.h>

//#include"Compute_stats.h"
//#include"draw_polygon_soup.h"

namespace PMP=CGAL::Polygon_mesh_processing;

////////////////////////////////////////////////////////////////////////////////
/// Build a polygon soup from an LCC.
/// op_bool = 1 -> Union
/// op_bool = 2 -> Intersection
/// op_bool = 3 -> A - B
/// op_bool = 4 -> B - A 
template<typename LCC>
void lcc_to_polygon_soup(const LCC& lcc,
                         std::vector<typename LCC::Point>& pts,
                         std::vector<std::vector<std::size_t>>& polygons)
{
  using VAD=typename LCC::Vertex_attribute_const_descriptor;
  std::map<VAD, std::size_t> lcc_to_pts;
  for(auto it=lcc.template one_dart_per_cell<2>().begin(),
        itend=lcc.template one_dart_per_cell<2>().end(); it!=itend; ++it)
  {
    typename LCC::Dart_const_descriptor dd= lcc.dart_descriptor(*it);
    polygons.push_back(std::vector<std::size_t>());
    std::vector<std::size_t>& p=polygons.back();
    do
    {
      auto res=lcc_to_pts.find(lcc.vertex_attribute(dd));
      if(res==lcc_to_pts.end())
      {
        res=lcc_to_pts.insert(std::make_pair(lcc.vertex_attribute(dd),
                                            pts.size())).first;
        auto point = lcc.point(dd);
        pts.push_back(point);
      }
      p.push_back(res->second);
      dd=lcc.template beta<1>(dd);
    }
    while(dd!=lcc.dart_descriptor(*it));
  }
}
////////////////////////////////////////////////////////////////////////////////
template<typename LCC>
class Polygon_soup_to_lcc_tool
{
  using FT=typename LCC::FT;
  using Point=typename LCC::Point;
  using Kernel=typename CGAL::Kernel_traits<Point>::Kernel;
  using VAD=typename LCC::Vertex_attribute_descriptor;
  using Dart_descriptor=typename LCC::Dart_descriptor;

public:
  Polygon_soup_to_lcc_tool(LCC& alcc,
                           const std::vector<Point>& pts,
                           const std::vector<std::vector<std::size_t>>& triangles)
      : m_lcc(alcc)
  {
    //std::cout<<"Create triangles and store edges..."<<std::flush;
    for (std::size_t i=0; i<triangles.size(); ++i)
    {
      VAD vad1=get_or_create_point(pts[triangles[i][0]]);
      VAD vad2=get_or_create_point(pts[triangles[i][1]]);
      VAD vad3=get_or_create_point(pts[triangles[i][2]]);

      Dart_descriptor dd1=m_lcc.make_triangle(vad1, vad2, vad3);
      Dart_descriptor dd2=m_lcc.make_triangle(vad2, vad1, vad3);
      m_lcc.template topo_sew<3>(dd1, dd2);

      // 3edges: (1) vad1->vad2 (2) vad2->vad3 (3) vad3->vad1
      // edges are stored using the smaller descriptor in first
      store_edge(dd1);
      store_edge(m_lcc.template beta<1>(dd1));
      store_edge(m_lcc.template beta<0>(dd1));
    }

    //std::cout<<std::endl<<"Sort edges..."<<std::flush;
    sort_edges();
    //std::cout<<std::endl<<"Sew2 edges..."<<std::flush;
    sew2_edges();
    //std::cout<<std::endl;

    // display_stats(m_lcc, true); // does not compile with epeck
  }

  VAD get_or_create_point(const Point& p)
  {
    auto it=m_ptov.find(p);
    if(it==m_ptov.end())
    { it=m_ptov.insert(std::make_pair(p, m_lcc.create_vertex_attribute(p))).first; }
    return it->second;
  }

  void store_edge(Dart_descriptor dd)
  {
    VAD vad1=m_lcc.vertex_attribute(dd);
    VAD vad2=m_lcc.vertex_attribute(m_lcc.template beta<1>(dd));

    if(vad1==vad2) { return; }
    if(vad1>vad2)
    {
      std::swap(vad1, vad2);
      dd=m_lcc.template beta<3>(dd);
    }

    auto it=m_edges.find(std::make_pair(vad1, vad2));
    if(it==m_edges.end())
    {
      it=m_edges.insert(std::make_pair(std::make_pair(vad1, vad2),
                                       std::vector<Dart_descriptor>())).first;
    }
    it->second.push_back(dd);
  }

  void sort_edges()
  {
    namespace pred=CGAL::Polygon_mesh_processing::Corefinement;

    for(auto& e:m_edges)
    {
      /* std::cout<<"Triangles around edge: ["<<m_lcc.point_of_vertex_attribute(e.first.first)
              <<" -> "<<m_lcc.point_of_vertex_attribute(e.first.second)<<"]: ";
      for(auto dd: e.second)
      { std::cout<<"("<<m_lcc.point(dd)<<"; "
                <<m_lcc.point(m_lcc.template beta<1>(dd))<<"; "
               <<m_lcc.point(m_lcc.template beta<0>(dd))<<")  "; }
      std::cout<<std::endl; */
     /* if(e.second.size()==1)
      {
        std::cout<<"[ERROR] in sort_edges: only one dart for an edge."<<std::endl;
      }
      else */if(e.second.size()>2)
      {
        std::vector<Dart_descriptor>& triangles=e.second;
        Point& a=m_lcc.point(triangles.front());
        Point& b=m_lcc.point(m_lcc.template beta<1>(triangles.front()));
        Point& c=m_lcc.point(m_lcc.template beta<0>(triangles.front()));
        auto less = [this, &a, &b, &c](Dart_descriptor d1, Dart_descriptor d2)
        {
          return pred::sorted_around_edge<Kernel>
              (a, b, c,
               m_lcc.point(m_lcc.template beta<0>(d1)),
               m_lcc.point(m_lcc.template beta<0>(d2)));
        };
        std::sort(triangles.begin()+1, triangles.end(), less);
      }
    }
  }

  void sew2_around_edge_dir1(std::vector<Dart_descriptor>& e)
  {
    Dart_descriptor dprev=e[0];
    Dart_descriptor dcur;
    std::size_t nb=1;
    for(; nb<e.size(); ++nb)
    {
      dcur=m_lcc.template beta<3>(e[nb]);

      assert(m_lcc.vertex_attribute(dprev)==
             m_lcc.vertex_attribute(m_lcc.template beta<1>(dcur)));
      assert(m_lcc.vertex_attribute(m_lcc.template beta<1>(dprev))==
             m_lcc.vertex_attribute(dcur));

      m_lcc.template topo_sew<2>(dprev, dcur);

      if(m_lcc.template is_free<2>(m_lcc.template beta<3,1>(dcur)))
      { m_totreat.push(m_lcc.template beta<3,1>(dcur)); }
      if(m_lcc.template is_free<2>(m_lcc.template beta<3,0>(dcur)))
      { m_totreat.push(m_lcc.template beta<3,0>(dcur)); }

      dprev=m_lcc.template beta<3>(dcur);
    }

    dcur=m_lcc.template beta<3>(e[0]);

    assert(m_lcc.vertex_attribute(dprev)==
           m_lcc.vertex_attribute(m_lcc.template beta<1>(dcur)));
    assert(m_lcc.vertex_attribute(m_lcc.template beta<1>(dprev))==
           m_lcc.vertex_attribute(dcur));

    m_lcc.template topo_sew<2>(dprev, dcur);

    if(m_lcc.template is_free<2>(m_lcc.template beta<3,1>(dcur)))
    { m_totreat.push(m_lcc.template beta<3,1>(dcur)); }
    if(m_lcc.template is_free<2>(m_lcc.template beta<3,0>(dcur)))
    { m_totreat.push(m_lcc.template beta<3,0>(dcur)); }
  }

  void sew2_around_edge_dir2(std::vector<Dart_descriptor>& e)
  {
    Dart_descriptor dprev=m_lcc.template beta<3>(e[0]);
    Dart_descriptor dcur;
    std::size_t nb=e.size()-1;
    for(; nb>0; --nb)
    {
      dcur=e[nb];

      assert(m_lcc.vertex_attribute(dprev)==
             m_lcc.vertex_attribute(m_lcc.template beta<1>(dcur)));
      assert(m_lcc.vertex_attribute(m_lcc.template beta<1>(dprev))==
             m_lcc.vertex_attribute(dcur));

      m_lcc.template topo_sew<2>(dprev, dcur);

      if(m_lcc.template is_free<2>(m_lcc.template beta<3,1>(dcur)))
      { m_totreat.push(m_lcc.template beta<3,1>(dcur)); }
      if(m_lcc.template is_free<2>(m_lcc.template beta<3,0>(dcur)))
      { m_totreat.push(m_lcc.template beta<3,0>(dcur)); }

      dprev=m_lcc.template beta<3>(dcur);
    }

    dcur=e[0];

    assert(m_lcc.vertex_attribute(dprev)==
           m_lcc.vertex_attribute(m_lcc.template beta<1>(dcur)));
    assert(m_lcc.vertex_attribute(m_lcc.template beta<1>(dprev))==
           m_lcc.vertex_attribute(dcur));

    m_lcc.template topo_sew<2>(dprev, dcur);

    if(m_lcc.template is_free<2>(m_lcc.template beta<3,1>(dcur)))
    { m_totreat.push(m_lcc.template beta<3,1>(dcur)); }
    if(m_lcc.template is_free<2>(m_lcc.template beta<3,0>(dcur)))
    { m_totreat.push(m_lcc.template beta<3,0>(dcur)); }
  }

  void sew2_around_edge_correct_orientation(Dart_descriptor dd)
  {
    if(!m_lcc.template is_free<2>(dd)) { return; }

    VAD va1=m_lcc.vertex_attribute(dd);
    VAD va2=m_lcc.vertex_attribute(m_lcc.template beta<1>(dd));
    if(va1<va2)
    {
      auto it=m_edges.find(std::make_pair(va1, va2));
      assert(it!=m_edges.end());
      sew2_around_edge_dir1(it->second);
    }
    else
    {
      auto it=m_edges.find(std::make_pair(va2, va1));
      assert(it!=m_edges.end());
      sew2_around_edge_dir2(it->second);
    }
  }

  void sew2_edges()
  {
    for(auto& e:m_edges)
    {
      // 1) Take an edge 2-free
      if(m_lcc.template is_free<2>(e.second[0]))
      {
        // 2) 2-sew all the triangles around the edge
        sew2_around_edge_dir1(e.second);

        // 3) iterate through adjacency triangles to ensure a consitent orientation
        while(!m_totreat.empty())
        {
          Dart_descriptor dd=m_totreat.front();
          m_totreat.pop();
          sew2_around_edge_correct_orientation(dd);
        }
      }
    }
  }

private:
  LCC& m_lcc;
  std::map<Point, VAD> m_ptov;
  std::map<std::pair<VAD, VAD>, std::vector<Dart_descriptor>> m_edges;
  std::queue<Dart_descriptor> m_totreat;
};


////////////////////////////////////////////////////////////////////////////////
/// Insert all the polygon in the LCC.
template<typename LCC>
void polygon_soup_to_lcc(std::vector<typename LCC::Point>& pts,
                         std::vector<std::vector<std::size_t>>& polygons,
                         LCC& lcc,
                         bool triangulate=true,
                         bool repair=true,
                         bool autorefine=true)
{
  if(triangulate)
  {
    std::cout<<"Triangulate polygons..."<<std::flush;
    PMP::triangulate_polygons(pts, polygons);
    std::cout<<"OK."<<std::endl;
    // draw(pts, polygons);
  }

  if(repair)
  {
    std::cout<<"Repair polygon soup..."<<std::flush;
    PMP::repair_polygon_soup(pts, polygons);
    std::cout<<"OK."<<std::endl;
    // draw(pts, polygons);
  }

  if(autorefine)
  {
    std::cout<<"Autorefine triangle soup..."<<std::flush;
    PMP::autorefine_triangle_soup(pts, polygons,
                                  CGAL::parameters::concurrency_tag
                                  (CGAL::Parallel_if_available_tag()));
    PMP::repair_polygon_soup(pts, polygons);
    std::cout<<"OK."<<std::endl;
    // draw(pts, polygons);
  }

  std::cout<<"Triangles to LCC..."<<std::flush;
  Polygon_soup_to_lcc_tool<LCC> psl(lcc, pts, polygons);
  std::cout<<"OK."<<std::endl;
}

#endif // LCC_TO_POLYGON_SOUP_H
////////////////////////////////////////////////////////////////////////////////
