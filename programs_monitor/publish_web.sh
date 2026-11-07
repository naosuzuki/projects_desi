#!/bin/bash
# Publish the monitor page to the FSU web server every minute (Nao 2026-11-04):
#   https://astrophysics.physics.fsu.edu/~suzuki/pca/   (the web host keeps its own home, not astro4's)
# monitor.html fetches status.json relatively, so both files go to ~/public_html/pca there.
#   nohup ./publish_web.sh > ~/data/sdssdesi_status/publish_web.log 2>&1 &
S=~/data/sdssdesi_status; D=public_html/pca
SSH="-o BatchMode=yes -o ConnectTimeout=20 -o PubkeyAcceptedAlgorithms=+ssh-rsa -o HostKeyAlgorithms=+ssh-rsa"
H=suzuki@astrophysics.physics.fsu.edu
ssh $SSH $H "mkdir -p $D && chmod 755 $D"
while true; do
  cp $S/monitor.html $S/index.html
  scp -q $SSH $S/status.json $S/monitor.html $S/index.html $H:$D/ || echo "$(date) scp failed"
  sleep 60
done
