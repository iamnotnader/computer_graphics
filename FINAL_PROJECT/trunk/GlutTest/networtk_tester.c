     #include <arpa/inet.h>
     #include <netinet/in.h>
     #include <stdio.h>
     #include <sys/types.h>
     #include <sys/socket.h>
     #include <unistd.h>
 	  #include <string.h>  
	  #include <stdlib.h>
	
     #define BUFLEN 512
     #define NPACK 10 
	
   void diep(char *s)
   {
      perror(s);
   }
   
   int main(int argc, char** argv)
   {
      struct sockaddr_in si_me, si_other;
      int s_in, s_out, i;
      socklen_t slen=sizeof(si_other);
      char buf[BUFLEN];
   
   	printf("MY_PORT THEIR_PORT THEIR_IP\n");
   	uint16_t MY_PORT = atoi(argv[2]);
   	uint16_t THEIR_PORT = atoi(argv[3]);
   	char* THEIR_IP = argv[4];
   
      if ((s_in=socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP))==-1)
         diep("socket");
   
      memset((char *) &si_me, 0, sizeof(si_me));
      si_me.sin_family = AF_INET;
      si_me.sin_port = htons(MY_PORT);
      si_me.sin_addr.s_addr = htonl(INADDR_ANY);
      if (bind(s_in, (struct sockaddr*)&si_me, sizeof(si_me))==-1)
         diep("bind"); 
   	  
      memset((char *) &si_other, 0, sizeof(si_other));
      si_other.sin_family = AF_INET;
      si_other.sin_port = htons(THEIR_PORT);
      if (inet_aton(THEIR_IP, &si_other.sin_addr)==0) {
         fprintf(stderr, "inet_aton() failed\n");
         exit(1);
      }
   
      while (1)
      {
         if (strcmp(argv[1], "-c") == 0)
         {
            for (i=0; i<NPACK; i++) 
            {
               printf("Sending packet %d\n", i);
               sprintf(buf, "This is packet %d\n", i);
            
               if (sendto(s_in, buf, BUFLEN, 0, (struct sockaddr*) &si_other, slen)==-1)
                  diep("sendto()");
            }
         	printf("%d\n",rand());
         	printf("waiting...\n");
            for (i=0; i<NPACK; i++) 
            {
               if (recvfrom(s_in, buf, BUFLEN, 0, (struct sockaddr*)&si_other, &slen)==-1)
                  diep("recvfrom()");
               printf("Received packet from %s:%d\nData: %s\n\n", 
                  inet_ntoa(si_other.sin_addr), ntohs(si_other.sin_port), buf);
            }
         }
         
         if (strcmp(argv[1], "-s") == 0)
         {
            for (i=0; i<NPACK; i++) 
            {
            	printf("waiting\n");
               if (recvfrom(s_in, buf, BUFLEN, 0, (struct sockaddr*)&si_other, &slen)==-1)
                  diep("recvfrom()");
               printf("Received packet from %s:%d\nData: %s\n\n", 
                  inet_ntoa(si_other.sin_addr), ntohs(si_other.sin_port), buf);
            }
            sleep(1);
         	printf("%d\n",rand());
            for (i=0; i<NPACK; i++) 
            {
               printf("Sending packet %d\n", i);
               sprintf(buf, "This is packet %d\n", i);
            
               if (sendto(s_in, buf, BUFLEN, 0, (struct sockaddr*)&si_other, slen)==-1)
                  diep("sendto()");
            }
         }
      }
   
   
      
         
      
   	 
      //close(s);
      return 0;
   }
