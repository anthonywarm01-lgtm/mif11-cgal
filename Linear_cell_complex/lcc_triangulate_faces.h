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
#ifndef LCC_TRIANGULATE_FACES_H
#define LCC_TRIANGULATE_FACES_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Linear_cell_complex_operations.h>
///////////////////////////////////////////////////////////////////////////////
template<typename Face_handle>
bool is_external(Face_handle fh)
{ return fh->info().is_external; }

////////////////////////////////////////////////////////////////////////////////
  /*!
   * \brief Get the number of edges of the triangular face given by the handle
   * \c fh that exist in the \c LinearCellComplex object.
   *
   * \tparam CDT_ \c CGAL::Constrained_Delaunay_triangulation_2<Traits,Tds,Itag>
   * typename
   *
   * \param fh a handle to a triangular face of the templated constrained
   * Delaunay triangulation data structure
   *
   * \return the number 0 <= n <= 3 of existing edges
   */
template<typename Face_handle>
int number_of_existing_edge(Face_handle fh)
{
  unsigned res=0;
  for(int i=0; i<3; ++i)
  { if(fh->info().exist_edge[i]) ++res; }
  return res;
}
////////////////////////////////////////////////////////////////////////////////
  /*!
   * \brief Get the index of the first free edge of the triangular face given by
   * the handle \c fh.
   *
   * An edge is free if the associated \c LinearCellComplex object has no
   * corresponding 1-cell.
   *
   * \tparam CDT_ \c CGAL::Constrained_Delaunay_triangulation_2<Traits,Tds,Itag>
   * typename
   *
   * \pre The number of existing edge is either 1 or 2
   *
   * \param fh a handle to a triangular face of the templated constrained
   * Delaunay triangulation data structure
   *
   * \return the index 0 <= i <= 3 of the first free edge
   */
template<typename Face_handle>
int get_first_free_edge(Face_handle fh)
{
  CGAL_assertion( number_of_existing_edge(fh)==1 || number_of_existing_edge(fh)==2 );
  for(int i=0; i<3; ++i)
  { if(!fh->info().exist_edge[i]) return i; }

  CGAL_assertion(false);
  return -1;
}
////////////////////////////////////////////////////////////////////////////////
  /*!
   * \brief Get the index of the a free edge of the triangular face given by
   * the handle \c fh.
   *
   * An edge is free if the associated \c LinearCellComplex object has no
   * corresponding 1-cell.
   *
   * \tparam CDT_ \c CGAL::Constrained_Delaunay_triangulation_2<Traits,Tds,Itag>
   * typename
   *
   * \pre The number of existing edge is either 1 or 2
   *
   * \param fh a handle to a triangular face of the templated constrained
   * Delaunay triangulation data structure
   *
   * \return the index 0 <= i <= 3 of a free edge
   */
template<typename Face_handle>
int get_free_edge(Face_handle fh)
{
  CGAL_assertion( number_of_existing_edge(fh)==2 );
  for(int i=0; i<3; ++i)
  { if(!fh->info().exist_edge[i]) return i; }

  CGAL_assertion(false);
  return -1;
}

template<typename Vertex_handle>
bool is_edge_insertable_case_one(Vertex_handle vh1, Vertex_handle vh2)
{ return vh1->info().degree_two && vh2->info().degree_two; }

///////////////////////////////////////////////////////////////////////////////
/// Triangulate the face containing dart d1.
template<typename LCC>
void constrained_delaunay_triangulation(LCC &lcc, typename LCC::Dart_descriptor d1)
{
  struct Vertex_info
  {
    Vertex_info(): dh(LCC::null_descriptor), degree_two(true)
    {}

    typename LCC::Dart_descriptor dh;
    bool degree_two;
  };

  struct Face_info
  {
    Face_info(): is_external(true), is_process(false)
    { exist_edge[0]=false; exist_edge[1]=false; exist_edge[2]=false; }

    bool exist_edge[3];
    bool is_external;
    bool is_process;
  };

  typedef CGAL::Projection_traits_3<CGAL::Exact_predicates_inexact_constructions_kernel> P_traits;
  typedef CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, P_traits> Vb;

  typedef CGAL::Triangulation_face_base_with_info_2<Face_info,P_traits> Fb1;

  typedef CGAL::Constrained_triangulation_face_base_2<P_traits, Fb1>    Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                   TDS;
  typedef CGAL::Exact_predicates_tag                                    Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits, TDS,
                                                     Itag>              CDT;

  if(lcc.template beta<1,1,1>(d1)==d1) { return; } // The face is already triangulated

  typename LCC::Vector normal=-CGAL::compute_normal_of_cell_2(lcc,d1);
  P_traits cdt_traits(normal);
  CDT cdt(cdt_traits);

  typename CDT::Vertex_handle vh1=nullptr, vh2=nullptr;
  std::vector<typename LCC::Dart_descriptor> toremove;
  typename LCC::size_type face=lcc.get_new_mark();
  lcc.template mark_cell<2,2>(d1, face);

  //inserting the constraints edge by edge
  typename LCC::template Dart_of_orbit_range<1>::iterator
    it(lcc.template darts_of_orbit<1>(d1).begin());
   for(typename LCC::template Dart_of_orbit_range<1>::iterator
         itend(lcc.template darts_of_orbit<1>(d1).end()); it!=itend; ++it)
   {
     if(!lcc.template is_free<2>(it) &&
        lcc.is_marked(lcc.template beta<2>(it), face))
     {
       if(it<lcc.template beta<2>(it))
       { toremove.push_back(it); }
     }
     else
     {
       vh1=cdt.insert(lcc.point(it)); // Does not insert if the point already exists
       vh1->info().dh=it; // Here we know that this dart does not belong to an inner edge
       vh1->info().degree_two=true;

       vh2=cdt.insert(lcc.point(lcc.other_extremity(it))); // Does not insert if the point already exists
       // Here we do not update dh and degree_two; this will be done later for the
       // dart having vh2 as point.

       cdt.insert_constraint(vh1, vh2);
   }
   }
   CGAL_assertion(cdt.is_valid());
   lcc.template unmark_cell<2,2>(d1, face);
   lcc.free_mark(face);

   std::queue<typename CDT::Face_handle> face_queue;
   typename CDT::Face_handle face_internal=nullptr;

   // Search a first internal face
   face_queue.push(cdt.infinite_vertex()->face());
   while(!face_queue.empty())
   {
     typename CDT::Face_handle fh=face_queue.front();
     face_queue.pop();
     if(!fh->info().is_process && face_internal==nullptr)
     {
       fh->info().is_process = true;
       for(int i=0; i<3; ++i)
       {
         if(!cdt.is_constrained(std::make_pair(fh, i)))
         { face_queue.push(fh->neighbor(i)); }
         else // assert(face_internal==nullptr);
         { face_internal=fh->neighbor(i); }
       }
     }
   }

   if(face_internal==nullptr)
   { return; } // No internal face => nothing to triangulate
   // Here only external faces are marked

   // Now we mark all internal faces
   face_queue.push(face_internal);
   while(!face_queue.empty())
   {
     typename CDT::Face_handle fh=face_queue.front();
     face_queue.pop();
     if(!fh->info().is_process)
     {
       fh->info().is_process =true;
       fh->info().is_external=false;
       for(int i = 0; i <3; ++i)
       {
         if(!cdt.is_constrained(std::make_pair(fh, i)))
         {
           if(!fh->neighbor(i)->info().is_process)
           { face_queue.push(fh->neighbor(i)); }
         }
         else
         { fh->info().exist_edge[i]=true; }
       }
     }
   }

   for(typename CDT::All_faces_iterator fit=cdt.all_faces_begin(),
         fitend=cdt.all_faces_end(); fit!=fitend; ++fit)
   {
     if(!fit->info().is_external)
     {
       int nbe=number_of_existing_edge(fit);
       if(nbe==1)
       {
         // triangle is (vc, vb, va)
         int index=get_first_free_edge(fit); // index is one edge which does not exist in the face
         typename CDT::Face_handle opposite_fh=fit->neighbor(index); // the opposite face along this edge
         const typename CDT::Vertex_handle va=fit->vertex(cdt. cw(index)); // source of edge index
         const typename CDT::Vertex_handle vb=fit->vertex(cdt.ccw(index)); // target of edge index
         // const typename CDT::Vertex_handle vc=fit->vertex(index);
         /* std::cout<<"Triangle "<<va->point()<<" "<<vb->point()<<" "<<vc->point()<<" one edge, insertable? "<<
                 (is_edge_insertable_case_one(va, vb)?"TRUE":"FALSE")<<std::endl; */

         if(is_edge_insertable_case_one(va, vb)) // if va and vb are degree 2, we can add the edge directly since there is no ambiguity
         {
           typename LCC::Dart_descriptor dd1=va->info().dh;
           typename LCC::Dart_descriptor dd2=vb->info().dh;
           // std::cout<<"I1: "<<lcc.point(dd1)<<" -> "<<lcc.point(dd2)<<std::endl;
           typename LCC::Dart_descriptor ndart=lcc.insert_cell_1_in_cell_2(dd1, dd2);

           va->info().dh=lcc.template beta<2>(ndart);
           fit->info().exist_edge[index]=true;
           opposite_fh->info().exist_edge[cdt.mirror_index(fit,index)]=true;

           va->info().degree_two=false;
           vb->info().degree_two=false;

           // TODO modify vb->info().dh and maybe we can remove the loop around vertices

           assert(number_of_existing_edge(fit)==2);
           face_queue.push(fit);
           assert(!opposite_fh->info().is_external);
           if(number_of_existing_edge(opposite_fh)==2)
           { face_queue.push(opposite_fh); }
         }
       }
       else if (nbe==2)
       { face_queue.push(fit); }
       // nbe==3 => nothing to do; nbe=0 too (will be considered later)
     }
   }

   // Remove internal edges (after the possible insertion just above)
   for(auto ittoremove: toremove)
   { lcc.template remove_cell<1>(ittoremove); }

   while(!face_queue.empty())
   {
     typename CDT::Face_handle fh=face_queue.front();
     face_queue.pop();
     CGAL_assertion(number_of_existing_edge(fh)>=2); // i.e. ==2 or ==3
     CGAL_assertion(!fh->info().is_external);

     if(number_of_existing_edge(fh)==2)
     {
       int index=get_free_edge(fh);
       typename CDT::Face_handle opposite_fh=fh->neighbor(index);

       CGAL_assertion( !fh->info().exist_edge[index] );
       CGAL_assertion( !opposite_fh->info().
                       exist_edge[cdt.mirror_index(fh,index)] );
       // triangle is (vc, vb, va)
       const typename CDT::Vertex_handle va=fh->vertex(cdt. cw(index));
       const typename CDT::Vertex_handle vb=fh->vertex(cdt.ccw(index));
       const typename CDT::Vertex_handle vc=fh->vertex(index);

       typename LCC::Dart_descriptor dd1=lcc.null_descriptor;
       for(typename LCC::template Dart_of_cell_range<0, 2>::iterator
             iti=lcc.template darts_of_cell<0, 2>(va->info().dh).begin();
           dd1==lcc.null_descriptor && iti.cont(); ++iti)
       {
         if(lcc.point(lcc.template beta<0>(iti))==vc->point())
         { dd1=iti; }
       }

       typename LCC::Dart_descriptor dd2=lcc.null_descriptor;
       for(typename LCC::template Dart_of_cell_range<0, 2>::iterator
             iti=lcc.template darts_of_cell<0, 2>(vb->info().dh).begin();
           dd2==lcc.null_descriptor && iti.cont(); ++iti)
       {
         if(lcc.point(lcc.template beta<1>(iti))==vc->point())
         { dd2=iti; }
       }

       //       assert(((lcc.beta<0,0>(dd1)==dd2) || lcc.beta<1,1>(dd1)==dd2));
       assert(dd1!=lcc.null_descriptor && dd2!=lcc.null_descriptor);
       //if(dd1!=nullptr && dd2!=nullptr) // TODO transform this if into an assert
       {
         typename LCC::Dart_descriptor ndart=lcc.insert_cell_1_in_cell_2(dd1, dd2);
         //std::cout<<"I2: "<<lcc.point(dd1)<<" -> "<<lcc.point(dd2)<<std::endl;

       va->info().dh=lcc.template beta<2>(ndart);

       fh->info().exist_edge[index]=true;
       opposite_fh->info().exist_edge[cdt.mirror_index(fh,index)]=true;

       if(!opposite_fh->info().is_external &&
          number_of_existing_edge(opposite_fh)==2)
       { face_queue.push(opposite_fh); }
       }
       /* else
       {
         std::cout<<"No dart found for edge ("<<va->point()<<") -> ("
                 <<vb->point()<<") for triangle with vc="<<vc->point()<<std::endl;
         return;
       } */
     }
   }
}
///////////////////////////////////////////////////////////////////////////////
/// Triangulate all marked faces (each face having at least one marked dart)
template<typename LCC>
void triangulate_marked_faces(LCC &lcc, typename LCC::size_type amark)
{
  // We are going to call constrained_delaunay_triangulation several time for
  // a same face when it has all its darts marked. But only the first call will
  // triangulate the face, the other ones will do nothing since the faces are
  // already triangulated. It is maybe faster to mark darts of the face in order
  // to avoid these successive calls, but not sure since we must unmark the
  // marked darts in another loop.
  for(auto it=lcc.darts().begin(); it!=lcc.darts().end(); ++it)
  {
    if(lcc.is_marked(it, amark))
    { constrained_delaunay_triangulation(lcc, it); }
  }
}
///////////////////////////////////////////////////////////////////////////////
/// Triangulate all faces of the lcc.
template<typename LCC>
void triangulate_all_faces(LCC &lcc)
{
  for(auto it=lcc.darts().begin(); it!=lcc.darts().end(); ++it)
  { constrained_delaunay_triangulation(lcc, it); }
}
///////////////////////////////////////////////////////////////////////////////
#endif // LCC_TRIANGULATE_FACES_H
