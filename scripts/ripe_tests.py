from warriorpy.ripe import systems as rs

if __name__ == "__main__":
    rp = rs.request_pool()
    tm = rs.thread_manager(10,
