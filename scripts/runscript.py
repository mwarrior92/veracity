import metric_validation as mv

m, fmt = mv.get_domain_matrix(1493596832, duration=60*60*24*2,
        mindomsize=1, country_set=["JP"])
