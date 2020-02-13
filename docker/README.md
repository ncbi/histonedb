# How to build and run docker

## Build image

1. Create an image

``` sudo docker image build -t histdbimage:1.0.0 . ```

2. List all images

``` sudo docker images ```

3. Remove image

``` sudo docker rmi <imageID> ```

4. Remove all images

``` sudo docker rmi $(sudo docker images -a -q) ```

## Run docker

1. Run docker image

``` sudo docker container run -it --name histonedb_docker -v /home/l_singh/docker_app_test/hist_docker_test:/var/www/histonedb --mount type=bind,src=/home/l_singh/docker_app_test/database,dst=/var/lib/mysql histdbimage:1.0.0 ```

``` sudo docker container run -it --name histonedb_docker -v /home/l_singh/docker_app_test/hist_docker_test:/var/www/histonedb --mount type=bind,src=/home/l_singh/docker_app_test/database,dst=/var/lib/mysql histdbimage:1.0.0 -without_setup```

2. List docker containers

``` sudo docker ps -a ```

3. Start docker container

``` sudo docker start <containerID> ```

or

``` sudo docker start -ia <containerID> ```

4. Connect to running container in interactive mode

``` sudo docker exec -it <containerID> bash ```

5. Remove container

``` sudo docker rm <containerID> ```

6. Remove all containers

``` sudo docker rm $(sudo docker ps -a -q) ```

7. Remove last container

``` sudo docker rm $(sudo docker ps -l -q) ```

8. Create docker container

``` sudo docker create -it --name histonedb_docker -v /home/l_singh/docker_app_test/hist_docker_test:/var/www/histonedb --mount type=bind,src=/home/l_singh/docker_app_test/database,dst=/var/lib/mysql histdbimage:1.0.0 -without_setup ```
