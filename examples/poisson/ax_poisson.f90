! Copyright (c) 2021-2025, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
module ax_poisson
  use ax_product
  use utils, only : neko_error
  use num_types, only : rp, i8
  use coefs, only : coef_t
  use space, only : space_t
  use mesh, only : mesh_t
  use math, only : addcol4, glsum
  implicit none
  private

  type, public, extends(ax_t) :: ax_poisson_t
   contains
     procedure, nopass :: compute => ax_poisson_compute
     procedure, pass(this) :: compute_vector => ax_poisson_compute_vector
  end type ax_poisson_t
  real(kind=rp), allocatable :: A_matrix(:,:)
  ! (i,j,val) format for the matrix
  ! NOTE: there may be multiple entries for the same i,j.
  ! in such a case they should be summed together (think of neighbor contributions)
  integer(kind=i8), allocatable :: A_rows(:)
  integer(kind=i8), allocatable :: A_cols(:)
  real(kind=rp), allocatable :: A_vals(:)
  real(kind=rp), allocatable :: u_vec(:), w_vec(:)
  

contains

  subroutine ax_poisson_compute(w, u, coef, msh, Xh)
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp) :: ur(Xh%lx, Xh%lx, Xh%lx)
    real(kind=rp) :: us(Xh%lx, Xh%lx, Xh%lx)
    real(kind=rp) :: ut(Xh%lx, Xh%lx, Xh%lx)
    real(kind=rp) :: wur(Xh%lx, Xh%lx, Xh%lx)
    real(kind=rp) :: wus(Xh%lx, Xh%lx, Xh%lx)
    real(kind=rp) :: wut(Xh%lx, Xh%lx, Xh%lx)
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l, num_dofs, irow, icol, idof, nnz

    ! @todo don't assume lx = ly = lz
    associate( D => Xh%dx, Dt => Xh%dxt, &
         G11 => coef%G11, G22 => coef%G22, G33 => coef%G33, &
         G12 => coef%G12, G13 => coef%G13, G23 => coef%G23, &
         n => msh%nelv, lx => Xh%lx)

      ! This is the matrix-based implementation of a matrix-free operator
      ! On the first call, it builds the matrix before computing the matvec
      ! On subsequent calls, it only computes the matvec using that matrix
      if (.not. allocated(A_matrix)) then
         write(*,*)
         write(*,*) '------------------------------'
         write(*,*) 'Size for matrices D, Dt, G11:'
         write(*,*) size(D)
         write(*,*) size(Dt)
         write(*,*) size(G11)
         ! write(*,*) G11
         write(*,*) 'Number of elements ', n
         write(*,*) 'Number of local DOFs ', lx
         write(*,*) 'DOF size ', coef%dof%size()
         ! write(*,*) coef%dof%dof
         write(*,*) 'Allocating matrix meow meow ^-^'
         ! true number of dofs (number of rows/cols in A)
         num_dofs = int(glsum(coef%mult, coef%dof%size()), i8)
         allocate(A_matrix(1,1))

         ! storing the matrix in (i,j,val) format needs one entry per dof contribution
         ! this means we will need elems
         nnz = n*lx*lx*lx*lx*3
         write(*,*) 'nnz'
         write(*,*) nnz
         allocate(A_vals(nnz))
         A_vals = 0.0_rp
         allocate(A_rows(nnz))
         allocate(A_cols(nnz))
         A_rows = 0_i8
         A_cols = 0_i8
         idof = 1 ! fortran indexing

         ! allocate the u and w vectors
         allocate(u_vec(num_dofs))
         allocate(w_vec(num_dofs))
      endif

      idof = 1

      ! Loop over mesh elements
      do e = 1, n
         write(*,*) 'Loop index ', e
         ! Compute the action of the derivative operator (D u)
         ! TODO: this assumes the mesh is structured
         !       so it's missing Gij with i .ne. j
         do k = 1, lx
            do j = 1, lx
               do i = 1, lx
                  tmp = 0.0_rp
                  do l = 1, lx
                     A_vals(idof) = G11(i,j,k,e) * D(i,l)
                     A_rows(idof) = coef%dof%dof(i,j,k,e)
                     A_cols(idof) = coef%dof%dof(l,j,k,e)
                     idof = idof + 1
                     ! A_matrix(coef%dof%dof(i,j,k,e),coef%dof%dof(l,j,k,e)) = D(i,l) ! TODO: this looks right
                     tmp = tmp + D(i,l) * u(l,j,k,e)
                  end do
                  ur(i,j,k) = G11(i,j,k,e) * tmp
                  
                  tmp = 0.0_rp
                  do l = 1, lx
                     A_vals(idof) = G22(i,j,k,e) * D(j,l)
                     A_rows(idof) = coef%dof%dof(i,j,k,e)
                     A_cols(idof) = coef%dof%dof(i,l,k,e)
                     idof = idof + 1
                     ! A_matrix(coef%dof%dof(i,j,k,e),coef%dof%dof(i,l,k,e)) = D(j,l) ! TODO: this looks right
                     tmp = tmp + D(j,l) * u(i,l,k,e)
                  end do
                  us(i,j,k) =  G22(i,j,k,e) * tmp

                  tmp = 0.0_rp
                  do l = 1, lx
                     A_vals(idof) = G33(i,j,k,e) * D(k,l)
                     A_rows(idof) = coef%dof%dof(i,j,k,e)
                     A_cols(idof) = coef%dof%dof(i,j,l,e)
                     idof = idof + 1
                     ! A_matrix(coef%dof%dof(i,j,k,e),coef%dof%dof(i,j,l,e)) = D(k,l) ! TODO: this looks right
                     tmp = tmp + D(k,l) * u(i,j,l,e)
                  end do
                  ut(i,j,k) =  G33(i,j,k,e) * tmp
               end do
            end do
         end do

         ! Compute the action of the derivative transpose operator (D^T (G (D u)))
         ! (D^T_xi u)
         do k = 1, lx
            do j = 1, lx
               do i = 1, lx
                  tmp = 0.0_rp
                  do l = 1, lx
                     tmp = tmp + Dt(i,l) * ur(l,j,k)
                  end do
                  w(i,j,k,e) = tmp
               end do
            end do
         end do

         ! (D^T_eta u)
         do k = 1, lx
            do j = 1, lx
               do i = 1, lx
                  tmp = 0.0_rp
                  do l = 1, lx
                     tmp = tmp + Dt(j,l) * us(i,l,k)
                  end do
                  w(i,j,k,e) = w(i,j,k,e) + tmp
               end do
            end do
         end do

         ! (D^T_gamma u)
         do k = 1, lx
            do j = 1, lx
               do i = 1, lx
                  tmp = 0.0_rp
                  do l = 1, lx
                     tmp = tmp + Dt(k,l) * ut(i,j,l)
                  end do
                  w(i,j,k,e) = w(i,j,k,e) + tmp
               end do
            end do
         end do

      end do ! e = 1, n

      write(*,*) 'idof ', idof

      u_vec = 0.0_rp
      w_vec = 0.0_rp

      ! Turn u and w into true vectors instead of the previous abominations
      do e = 1, n
         do i = 1, lx
            do j = 1, lx
               do k = 1, lx
                  u_vec(coef%dof%dof(i,j,k,e)) = u_vec(coef%dof%dof(i,j,k,e)) + u(i,j,k,e) * coef%mult(i,j,k,e)
                  w_vec(coef%dof%dof(i,j,k,e)) = w_vec(coef%dof%dof(i,j,k,e)) + w(i,j,k,e) * coef%mult(i,j,k,e)
               end do
            end do
         end do
      end do

      ! Do the true matrix-vector multiplication, (i,j,val) style
      do i = 1, nnz
         irow = A_rows(i)
         icol = A_cols(i)
         w_vec(irow) = w_vec(irow) + A_vals(i) * u_vec(icol)
      end do

    end associate
  end subroutine ax_poisson_compute

  subroutine ax_poisson_compute_vector(this, au, av, aw, u, v, w, coef, msh, Xh)
    class(ax_poisson_t), intent(in) :: this
    type(space_t), intent(in) :: Xh
    type(mesh_t), intent(in) :: msh
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(inout) :: au(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: av(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: aw(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: v(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)

    call neko_error('Not in use')

  end subroutine ax_poisson_compute_vector

end module ax_poisson
