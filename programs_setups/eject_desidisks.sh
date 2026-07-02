#!/usr/bin/env bash
#
# eject_desidisks.sh
# ------------------
# One-shot clean pause + eject for the DESI download disks.
# Stops Globus Connect Personal (pauses the transfer, no progress lost) — which
# is the only thing holding the disks open — then ejects the external
# "My Book" drives so they're safe to unplug.
#
# The Globus task resumes automatically next time Globus Connect Personal runs.
#
# NOTE: deliberately does NOT use `lsof +D` — that recursively scans the whole
# multi-TB disk (slow) and the scan itself blocks the unmount. Stopping Globus
# Connect Personal releases the disk; diskutil eject then just works.
#
# Usage:  ./eject_desidisks.sh
#
set -uo pipefail

echo "1) stopping Globus Connect Personal (pauses transfer)..."
pkill -9 -f "Globus Connect Personal" 2>/dev/null
pkill -9 -f "globus-gridftp-server"   2>/dev/null
pkill -9 -f "relaytool"               2>/dev/null
sleep 3

echo "2) ejecting external disks..."
ok=1
for v in "/Volumes/My Book 1" "/Volumes/My Book"; do
    [ -d "$v" ] || continue                      # skip if not mounted
    if diskutil eject "$v" >/dev/null 2>&1; then
        echo "   EJECTED: $v"
    elif diskutil unmount force "$v" >/dev/null 2>&1; then
        echo "   FORCE-UNMOUNTED: $v"
    else
        echo "   FAILED: $v (something still has it open)"; ok=0
    fi
done

echo "3) result:"
if df -h 2>/dev/null | grep -qiE "/Volumes/My Book"; then
    df -h | grep -iE "/Volumes/My Book" | sed 's/^/   still mounted: /'
    ok=0
else
    echo "   both external disks unmounted — safe to unplug"
fi

[ "$ok" = 1 ] && echo "Done." || { echo "Re-run, or check the holder with: diskutil eject '/Volumes/My Book 1'"; exit 1; }
