!
!  MPI alias definitions for serial machines.
!
!
!  last revised: 15 January 2011
!
      module mpidefs
      implicit none
      integer :: mpi_comm_world,ms_mpi_comm_world,ms_mpi_sum,ms_mpi_max,ms_mpi_min

      contains

         subroutine ms_mpi(mpi_command,mpi_recv_buf_i,mpi_recv_buf_r,mpi_recv_buf_c,mpi_recv_buf_dp, &
                           mpi_recv_buf_dc,mpi_send_buf_i,mpi_send_buf_r,mpi_send_buf_c, &
                           mpi_send_buf_dp,mpi_send_buf_dc,mpi_number,mpi_comm,mpi_group,mpi_rank,mpi_size, &
                           mpi_new_comm,mpi_new_group,mpi_new_group_list,mpi_operation)
         integer, optional :: mpi_number,mpi_recv_buf_i(*),mpi_send_buf_i(*),mpi_comm,mpi_group,mpi_rank, &
                              mpi_size,mpi_new_comm,mpi_new_group,mpi_new_group_list(*),mpi_operation
         integer :: stat(1)
         real(4), optional :: mpi_recv_buf_r(*),mpi_send_buf_r(*)
         real(8), optional :: mpi_recv_buf_dp(*),mpi_send_buf_dp(*)
         complex(4), optional :: mpi_recv_buf_c(*),mpi_send_buf_c(*)
         complex(8), optional :: mpi_recv_buf_dc(*),mpi_send_buf_dc(*)
         character(*) :: mpi_command
         integer :: type,ierr,comm,size,rank,group,newcomm

         if(mpi_command.eq.'init') then
            mpi_comm_world=1
            return
         endif
         if(mpi_command.eq.'finalize') then
            return
         endif
         if(present(mpi_comm)) then
            comm=mpi_comm
         else
            comm=mpi_comm_world
         endif
         if(mpi_command.eq.'size') then
            mpi_size=1
            return
         endif
         if(mpi_command.eq.'rank') then
            mpi_rank=0
            return
         endif
         if(mpi_command.eq.'group') then
            return
         endif
         if(mpi_command.eq.'incl') then
            mpi_new_group=0
            return
         endif
         if(mpi_command.eq.'create') then
            mpi_new_comm=1
            return
         endif
         if(mpi_command.eq.'barrier') then
            return
         endif

         if(present(mpi_recv_buf_i).or.present(mpi_send_buf_i)) then
            if(mpi_command.eq.'bcast'.or.mpi_command.eq.'send'.or.mpi_command.eq.'recv') then
               return
            endif
            if(mpi_command.eq.'reduce'.or.mpi_command.eq.'allreduce') then
               if(present(mpi_send_buf_i)) then
                  mpi_recv_buf_i(1:mpi_number)=mpi_send_buf_i(1:mpi_number)
               endif
            endif
            return
         endif

         if(present(mpi_recv_buf_r).or.present(mpi_send_buf_r)) then
            if(mpi_command.eq.'bcast'.or.mpi_command.eq.'send'.or.mpi_command.eq.'recv') then
               return
            endif
            if(mpi_command.eq.'reduce'.or.mpi_command.eq.'allreduce') then
               if(present(mpi_send_buf_r)) then
                  mpi_recv_buf_r(1:mpi_number)=mpi_send_buf_r(1:mpi_number)
               endif
            endif
            return
         endif

         if(present(mpi_recv_buf_c).or.present(mpi_send_buf_c)) then
            if(mpi_command.eq.'bcast'.or.mpi_command.eq.'send'.or.mpi_command.eq.'recv') then
               return
            endif
            if(mpi_command.eq.'reduce'.or.mpi_command.eq.'allreduce') then
               if(present(mpi_send_buf_c)) then
                  mpi_recv_buf_c(1:mpi_number)=mpi_send_buf_c(1:mpi_number)
               endif
            endif
            return
         endif

         if(present(mpi_recv_buf_dp).or.present(mpi_send_buf_dp)) then
            if(mpi_command.eq.'bcast'.or.mpi_command.eq.'send'.or.mpi_command.eq.'recv') then
               return
            endif
            if(mpi_command.eq.'reduce'.or.mpi_command.eq.'allreduce') then
               if(present(mpi_send_buf_dp)) then
                  mpi_recv_buf_dp(1:mpi_number)=mpi_send_buf_dp(1:mpi_number)
               endif
            endif
            return
         endif

         if(present(mpi_recv_buf_dc).or.present(mpi_send_buf_dc)) then
            if(mpi_command.eq.'bcast'.or.mpi_command.eq.'send'.or.mpi_command.eq.'recv') then
               return
            endif
            if(mpi_command.eq.'reduce'.or.mpi_command.eq.'allreduce') then
               if(present(mpi_send_buf_dc)) then
                  mpi_recv_buf_dc(1:mpi_number)=mpi_send_buf_dc(1:mpi_number)
               endif
            endif
            return
         endif

         end subroutine ms_mpi
      end module mpidefs
