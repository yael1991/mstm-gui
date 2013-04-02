!
!  MPI alias definitions for parallel machines.
!
!
!  last revised: 15 January 2011
!
      module mpidefs
      use mpi
      implicit none
      integer :: ms_mpi_comm_world,ms_mpi_sum,ms_mpi_max,ms_mpi_min

      contains

         subroutine ms_mpi(mpi_command,mpi_recv_buf_i,mpi_recv_buf_r,mpi_recv_buf_c,mpi_recv_buf_dp, &
                           mpi_recv_buf_dc,mpi_send_buf_i,mpi_send_buf_r,mpi_send_buf_c, &
                           mpi_send_buf_dp,mpi_send_buf_dc,mpi_number,mpi_comm,mpi_group,mpi_rank,mpi_size,&
                           mpi_new_comm,mpi_new_group,mpi_new_group_list,mpi_operation)
         integer, optional :: mpi_number,mpi_recv_buf_i(*),mpi_send_buf_i(*),mpi_comm,mpi_group,mpi_rank, &
                              mpi_size,mpi_new_comm,mpi_new_group,mpi_new_group_list(*),mpi_operation
         integer :: stat(MPI_STATUS_SIZE)
         real(4), optional :: mpi_recv_buf_r(*),mpi_send_buf_r(*)
         real(8), optional :: mpi_recv_buf_dp(*),mpi_send_buf_dp(*)
         complex(4), optional :: mpi_recv_buf_c(*),mpi_send_buf_c(*)
         complex(8), optional :: mpi_recv_buf_dc(*),mpi_send_buf_dc(*)
         character(*) :: mpi_command
         integer :: type,ierr,comm,size,rank,group,newcomm

         if(mpi_command.eq.'init') then
            call mpi_init(ierr)
            ms_mpi_comm_world=mpi_comm_world
            ms_mpi_sum=mpi_sum
            ms_mpi_max=mpi_max
            ms_mpi_min=mpi_min
            return
         endif
         if(mpi_command.eq.'finalize') then
            call mpi_finalize(ierr)
            return
         endif
         if(present(mpi_comm)) then
            comm=mpi_comm
         else
            comm=mpi_comm_world
         endif
         if(mpi_command.eq.'size') then
            call mpi_comm_size(comm,size,ierr)
            mpi_size=size
            return
         endif
         if(mpi_command.eq.'rank') then
            call mpi_comm_rank(comm,rank,ierr)
            mpi_rank=rank
            return
         endif
         if(mpi_command.eq.'group') then
            call mpi_comm_group(comm,group,ierr)
            mpi_group=group
            return
         endif
         if(mpi_command.eq.'incl') then
            call mpi_group_incl(mpi_group,mpi_size,mpi_new_group_list,group,ierr)
            mpi_new_group=group
            return
         endif
         if(mpi_command.eq.'create') then
            call mpi_comm_create(comm,mpi_group,newcomm,ierr)
            mpi_new_comm=newcomm
            return
         endif
         if(mpi_command.eq.'barrier') then
            call mpi_barrier (comm,ierr)
            return
         endif

         if(present(mpi_recv_buf_i).or.present(mpi_send_buf_i)) then
            type=mpi_integer
            if(mpi_command.eq.'bcast') then
               call MPI_BCAST (mpi_send_buf_i,mpi_number,type,mpi_rank,comm,ierr)
               return
            endif
            if(mpi_command.eq.'send') then
               call mpi_send(mpi_send_buf_i,mpi_number,type,mpi_rank,1,comm,ierr)
               return
            endif
            if(mpi_command.eq.'recv') then
               call mpi_recv(mpi_recv_buf_i,mpi_number,type,mpi_rank,1,comm,stat,ierr)
               return
            endif
            if(mpi_command.eq.'reduce') then
               if(present(mpi_send_buf_i)) then
                  call mpi_reduce(mpi_send_buf_i,mpi_recv_buf_i,mpi_number,type,mpi_operation, &
                               mpi_rank,comm,ierr)
               else
                  call mpi_reduce(mpi_in_place,mpi_recv_buf_i,mpi_number,type,mpi_operation, &
                               mpi_rank,comm,ierr)
               endif
               return
            endif
            if(mpi_command.eq.'allreduce') then
               if(present(mpi_send_buf_i)) then
                  call mpi_allreduce(mpi_send_buf_i,mpi_recv_buf_i,mpi_number,type,mpi_operation, &
                               comm,ierr)
               else
                  call mpi_allreduce(mpi_in_place,mpi_recv_buf_i,mpi_number,type,mpi_operation, &
                               comm,ierr)
               endif
               return
            endif
         endif

         if(present(mpi_recv_buf_r).or.present(mpi_send_buf_r)) then
            type=mpi_real
            if(mpi_command.eq.'bcast') then
               call MPI_BCAST (mpi_send_buf_r,mpi_number,type,mpi_rank,comm,ierr)
               return
            endif
            if(mpi_command.eq.'send') then
               call mpi_send(mpi_send_buf_r,mpi_number,type,mpi_rank,1,comm,ierr)
               return
            endif
            if(mpi_command.eq.'recv') then
               call mpi_recv(mpi_recv_buf_r,mpi_number,type,mpi_rank,1,comm,stat,ierr)
               return
            endif
            if(mpi_command.eq.'reduce') then
               if(present(mpi_send_buf_r)) then
                  call mpi_reduce(mpi_send_buf_r,mpi_recv_buf_r,mpi_number,type,mpi_operation, &
                               mpi_rank,comm,ierr)
               else
                  call mpi_reduce(mpi_in_place,mpi_recv_buf_r,mpi_number,type,mpi_operation, &
                               mpi_rank,comm,ierr)
               endif
               return
            endif
            if(mpi_command.eq.'allreduce') then
               if(present(mpi_send_buf_r)) then
                  call mpi_allreduce(mpi_send_buf_r,mpi_recv_buf_r,mpi_number,type,mpi_operation, &
                               comm,ierr)
               else
                  call mpi_allreduce(mpi_in_place,mpi_recv_buf_r,mpi_number,type,mpi_operation, &
                               comm,ierr)
               endif
               return
            endif
         endif

         if(present(mpi_recv_buf_c).or.present(mpi_send_buf_c)) then
            type=mpi_complex
            if(mpi_command.eq.'bcast') then
               call MPI_BCAST (mpi_send_buf_c,mpi_number,type,mpi_rank,comm,ierr)
               return
            endif
            if(mpi_command.eq.'send') then
               call mpi_send(mpi_send_buf_c,mpi_number,type,mpi_rank,1,comm,ierr)
               return
            endif
            if(mpi_command.eq.'recv') then
               call mpi_recv(mpi_recv_buf_c,mpi_number,type,mpi_rank,1,comm,stat,ierr)
               return
            endif
            if(mpi_command.eq.'reduce') then
               if(present(mpi_send_buf_c)) then
                  call mpi_reduce(mpi_send_buf_c,mpi_recv_buf_c,mpi_number,type,mpi_operation, &
                               mpi_rank,comm,ierr)
               else
                  call mpi_reduce(mpi_in_place,mpi_recv_buf_c,mpi_number,type,mpi_operation, &
                               mpi_rank,comm,ierr)
               endif
               return
            endif
            if(mpi_command.eq.'allreduce') then
               if(present(mpi_send_buf_c)) then
                  call mpi_allreduce(mpi_send_buf_c,mpi_recv_buf_c,mpi_number,type,mpi_operation, &
                               comm,ierr)
               else
                  call mpi_allreduce(mpi_in_place,mpi_recv_buf_c,mpi_number,type,mpi_operation, &
                               comm,ierr)
               endif
               return
            endif
         endif

         if(present(mpi_recv_buf_dp).or.present(mpi_send_buf_dp)) then
            type=mpi_double_precision
            if(mpi_command.eq.'bcast') then
               call MPI_BCAST (mpi_send_buf_dp,mpi_number,type,mpi_rank,comm,ierr)
               return
            endif
            if(mpi_command.eq.'send') then
               call mpi_send(mpi_send_buf_dp,mpi_number,type,mpi_rank,1,comm,ierr)
               return
            endif
            if(mpi_command.eq.'recv') then
               call mpi_recv(mpi_recv_buf_dp,mpi_number,type,mpi_rank,1,comm,stat,ierr)
               return
            endif
            if(mpi_command.eq.'reduce') then
               if(present(mpi_send_buf_dp)) then
                  call mpi_reduce(mpi_send_buf_dp,mpi_recv_buf_dp,mpi_number,type,mpi_operation, &
                               mpi_rank,comm,ierr)
               else
                  call mpi_reduce(mpi_in_place,mpi_recv_buf_dp,mpi_number,type,mpi_operation, &
                               mpi_rank,comm,ierr)
               endif
               return
            endif
            if(mpi_command.eq.'allreduce') then
               if(present(mpi_send_buf_dp)) then
                  call mpi_allreduce(mpi_send_buf_dp,mpi_recv_buf_dp,mpi_number,type,mpi_operation, &
                               comm,ierr)
               else
                  call mpi_allreduce(mpi_in_place,mpi_recv_buf_dp,mpi_number,type,mpi_operation, &
                               comm,ierr)
               endif
               return
            endif
         endif

         if(present(mpi_recv_buf_dc).or.present(mpi_send_buf_dc)) then
            type=mpi_double_complex
            if(mpi_command.eq.'bcast') then
               call MPI_BCAST (mpi_send_buf_dc,mpi_number,type,mpi_rank,comm,ierr)
               return
            endif
            if(mpi_command.eq.'send') then
               call mpi_send(mpi_send_buf_dc,mpi_number,type,mpi_rank,1,comm,ierr)
               return
            endif
            if(mpi_command.eq.'recv') then
               call mpi_recv(mpi_recv_buf_dc,mpi_number,type,mpi_rank,1,comm,stat,ierr)
               return
            endif
            if(mpi_command.eq.'reduce') then
               if(present(mpi_send_buf_dc)) then
                  call mpi_reduce(mpi_send_buf_dc,mpi_recv_buf_dc,mpi_number,type,mpi_operation, &
                               mpi_rank,comm,ierr)
               else
                  call mpi_reduce(mpi_in_place,mpi_recv_buf_dc,mpi_number,type,mpi_operation, &
                               mpi_rank,comm,ierr)
               endif
               return
            endif
            if(mpi_command.eq.'allreduce') then
               if(present(mpi_send_buf_dc)) then
                  call mpi_allreduce(mpi_send_buf_dc,mpi_recv_buf_dc,mpi_number,type,mpi_operation, &
                               comm,ierr)
               else
                  call mpi_allreduce(mpi_in_place,mpi_recv_buf_dc,mpi_number,type,mpi_operation, &
                               comm,ierr)
               endif
               return
            endif
         endif

         end subroutine ms_mpi
      end module mpidefs
